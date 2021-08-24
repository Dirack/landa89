/* Landa 1988 experiment: VFSA velocity inversion based on stereotomography and NIP tomography strategies

This program is a reproduction of the experiment in the article 'A method for determination of velocity
and depth from seismic reflection data' avaliable in the doc directory of this repository

The initial velocity model and NIP sources position used in this program is set up using sfnipmodsetup. This program does the forward modelling by ray tracying from NIP sources to surface and gets reflection traveltime.

The time misfit is calculated by the difference between the reflection traveltime obtained in the forward modelling and the traveltime calculated by CRE traveltime approximation formula for RNIP and BETA parameters given. This time misfit is used as a convergence criteria for VFSA global optimization algorithm to obtain optimized velocity model.

*/

#include <math.h>
#include <rsf.h>
#include "velocity_lib.h"

int main(int argc, char* argv[])
{
	bool verb; // Verbose parameter
	int n[2]; // Velocity grid dimensions n[0]=n1, n[1]=n2
	float d[2]; // Velocity grid sampling d[0]=d1, d[1]=d2
	float o[2]; // Velocity grid origin o[0]=o1, o[1]=o2
	float** s; // NIP source position (z,x)
	float* cnewv; // Temporary parameters vector used in VFSA
	float* cnewz; // Temporary parameters vector used in VFSA
	float* otsv; // Optimized parameters vector
	float* otsz; // Optimized parameters vector
	float tmis0=100; // Best time misfit
	float otmis=0; // Best time misfit
	float deltaE; // Delta (Metrópolis criteria in VFSA)
	float Em0=0; // Energy (VFSA algorithm)
	float PM; // Metrópolis criteria
	float temp=1; // Temperature for VFSA algorithm
	float u=0; // Random number between 0 and 1
	int nit; // Number of VFSA iterations
	float temp0; // Initial temperature for VFSA
	float c0; // Damping factor for VFSA
	int ndim; // n1 dimension in shotsfile, should be equal 2
	int nshot; // n2 dimensions in shotsfile, number of shots
	int nm; // Number of samples in velocity grid n1*n2
	float* a; // Normal Ray initial angle for each NIP source
	float* slow; // slowness model
	int im; // loop counter
	float v; // Velocity temporary variable
	float v0; // Near surface velocity
	int ns; // Number of NIP sources
	int q; // Loop counter for VFSA iteration
	float tmis; // data time misfit value
	float *m0; // CMP's for normal rays
	float *t0; // t0's for normal rays
	float *RNIP; // Rnip parameters vector
	float *BETA; // Beta parameters vector
	float* sz; // Depth coordinates of the spline velocity function
	int nsz; // Dimension of sz vector
	float dsz;
	float osz;
	int nsv;
	float* sv; // Velocity coordinates of the spline velocity function
	sf_file shots; // NIP sources (z,x)
	sf_file vel; // background velocity model
	sf_file velinv; // Inverted velocity model
	sf_file angles; // Normal ray angles (degrees)
	sf_file m0s; // Central CMPs m0
	sf_file t0s; // Normal ray traveltimes
	sf_file rnips; // RNIP parameter for each m0
	sf_file betas; // BETA parameter for each m0
	sf_file sz_file; // z coordinates of the cubic spline functions
	sf_file vz_file; // v coordinates of the cubic spline functions
	sf_file vspline; // Cubic spline velocity model
	sf_file zspline;

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
	
	if(!sf_getbool("verb",&verb)) verb=true;
	/* verbose parameter (y/n) */

	if(!sf_getfloat("v0",&v0)) v0=1.5;
	/* Near surface velocity (Km/s) */

	if(!sf_getint("nit",&nit)) nit=1;
	/* Number of VFSA iterations */

	if(!sf_getfloat("temp0",&temp0)) temp0=5;
	/* Initial temperature for VFSA algorithm */

	if(!sf_getfloat("c0",&c0)) c0=0.1;
	/* Damping factor for VFSA algorithm */

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

	/*if(verb){
		sf_warning("Command line Parameters");
		sf_warning("v0=%f nit=%d temp0=%f c0=%f",v0,nit,temp0,c0);
		sf_warning("Input file (Velocity model)");
		sf_warning("n1=%d d1=%f o1=%f",*n,*d,*o);
		sf_warning("n2=%d d2=%f o2=%f",*(n+1),*(d+1),*(o+1));
		sf_warning("Input file (shotsfile)");
		sf_warning("n1=%d",ndim);
		sf_warning("n2=%d",nshot);
		sf_warning("Input file (anglefile)");
		sf_warning("n1=%d",ns);
		sf_warning("Input file (sz)");
		sf_warning("nz=%d",nsz);
	}*/

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

	/* Generate optimal velocity model */
	updateVelocityModel(n,o,d,sv,nsv,sz,nsz,osz,dsz,slow,nm);

	/* Write velocity model file */
	sf_floatwrite(slow,nm,velinv);
}
