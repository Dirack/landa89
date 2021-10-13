/* Build velocity model with a lateral velocity variation layer for landa89 experiment

TODO
*/

#include <math.h>
#include <rsf.h>
#include "velocity_lib.h"

int main(int argc, char* argv[])
{
	int n[2]; // Velocity grid dimensions n[0]=n1, n[1]=n2
	float d[2]; // Velocity grid sampling d[0]=d1, d[1]=d2
	float o[2]; // Velocity grid origin o[0]=o1, o[1]=o2
	int nm; // Number of samples in velocity grid n1*n2
	float* slow; // slowness model
	int im; // loop counter
	float v; // Velocity temporary variable
	float* sz; // Depth coordinates of the spline velocity function
	int nsz; // Dimension of sz vector
	float dsz;
	float osz;
	int nsv;
	float* sv; // Velocity coordinates of the spline velocity function
	sf_file vel; // background velocity model
	sf_file velinv; // Inverted velocity model
	sf_file sz_file; // z coordinates of the cubic spline functions
	sf_file vz_file; // v coordinates of the cubic spline functions

	sf_init(argc,argv);

	vel = sf_input("in");
	velinv = sf_output("out");
	sz_file = sf_input("sz");
	vz_file = sf_input("vz");

	/* Velocity model: get 2D grid parameters */
	if(!sf_histint(vel,"n1",n)) sf_error("No n1= in input");
	if(!sf_histint(vel,"n2",n+1)) sf_error("No n2= in input");
	if(!sf_histfloat(vel,"d1",d)) sf_error("No d1= in input");
	if(!sf_histfloat(vel,"d2",d+1)) sf_error("No d2= in input");
	if(!sf_histfloat(vel,"o1",o)) o[0]=0.;
	if(!sf_histfloat(vel,"o2",o+1)) o[1]=0.;
	
	/* Cubic spline vector */
	if(!sf_histint(sz_file,"n1",&nsz)) sf_error("No n1= in sz file");
	if(!sf_histfloat(sz_file,"d1",&dsz)) sf_error("No d1= in sz file");
	if(!sf_histfloat(sz_file,"o1",&osz)) sf_error("No o1= in sz file");
	if(!sf_histint(vz_file,"n1",&nsv)) sf_error("No n1= in sv file");

	/* Build cubic spline velocity matrix */
	sv = sf_floatalloc(nsv);
	sz = sf_floatalloc(nsz);
	sf_floatread(sz,nsz,sz_file);
	sf_floatread(sv,nsv,vz_file);

	/* get slowness squared (Background model) */
	nm = n[0]*n[1];
	slow =  sf_floatalloc(nm);
	sf_floatread(slow,nm,vel);

	for(im=0;im<nm;im++){
		v = slow[im];
		slow[im] = 1./(v*v);
	}

	/* Generate optimal velocity model */
	updateVelocityModel(n,o,d,sv,nsv,sz,nsz,osz,dsz,slow,nm);

	/* Velocity model from inversion */
	sf_putint(velinv,"n1",n[0]);
	sf_putint(velinv,"n2",n[1]);
	sf_putint(velinv,"n3",1);
	sf_putfloat(velinv,"d1",d[0]);
	sf_putfloat(velinv,"d2",d[1]);
	sf_putfloat(velinv,"o1",o[0]);
	sf_putfloat(velinv,"o2",o[1]);
	sf_putfloat(velinv,"d3",1);
	sf_putfloat(velinv,"o3",0);

	/* Write velocity model file */
	sf_floatwrite(slow,nm,velinv);
}
