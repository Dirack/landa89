/*
	 velocity_lib.c (c)
	 
	 Purpose: Functions to update velocity model.
	 	 
	 Version 1.0
	 
	 Site: https://dirack.github.io
	 
	 Programmer: Rodolfo A. C. Neves (Dirack) 15/10/2021

	 Email:  rodolfo_profissional@hotmail.com

	 License: GPL-3.0 <https://www.gnu.org/licenses/gpl-3.0.txt>.

*/

#include <stdio.h>
#include <stdlib.h>
#include <rsf.h>
#include "velocity_lib.h"

void calculateSplineCoeficients(int n, /* Vectors (x,y) dimension */
				float* x, /* x coordinates */
				float* y, /* y coordinates */
				float** coef, /* Spline coefficients */
				int n_stripes /* Model x axis is divided in n stripes */)
/*< Function to calculate natural cubic spline coefficients

Note: It Receives n points and two vectors x and y with n dimension.
It returns a coefficients vector with 4 coefficients for each of the
n-1 natural cubic splines, coef[(n-1)*4].

IMPORTANT: The number of points must be equal or major than 3 (n>3)
and x vector must be in crescent order.

>*/
{

	float s2[n]; // Second derivatives matrix
	int i, ip1, ip2, im1, m, k; // Loop counter
	float hb, ha, deltaa, deltab, t; // temporary variables
	float e[n-2]; // hi's vector
	float dp[n-2]; // main diagonal

	/* Vectors dimension must be major than 3 */
	if(n<3){
		fprintf(stderr,"Erro, n<3\n");
		exit(-1);
	}

	/* x vector must be in crescent order */
	for(i=1;i<n;i++){
		if(x[i-1]>x[i]){
			fprintf(stderr,"Erro, vetor x deve possuir ordem crescente\n");
			exit(-2);
		}
	}

	for(k=0;k<n_stripes;k++){
		
		/* Simetric tridiagonal linear system build */
		ha = x[1]-x[0]; deltaa = (y[k*n+1]-y[k*n+0])/ha; m=n-2;
		for(i=0;i<m;i++){
			ip1 = i+1; ip2 = i+2;
			hb = x[ip2]-x[ip1];
			deltab = (y[k*n+ip2]-y[k*n+ip1])/hb;
			e[i] = hb; dp[i] = 2*(ha+hb);
			s2[ip1] = 6*(deltab-deltaa);
			ha=hb; deltaa=deltab;
		}

		/* Gauss elimination */
		for(i=1;i<m;i++){
			ip1=i+1; im1=i-1;
			t = e[im1]/dp[im1];
			dp[i] = dp[i]-t*e[im1];
			s2[ip1] = s2[ip1]-t*s2[i];
		}

		/* Retroactive substitutive solution */
		s2[m]=s2[m]/dp[m-1];
		for(i=m-1;i>0;i--){
			ip1=i+1; im1=i-1;
			s2[i]=(s2[i]-e[im1]*s2[ip1])/dp[im1];
		}
		s2[0]=0; s2[n-1]=0;

		/* Calculate spline coefficients */
		for(i=0;i<n-1;i++){
			ha = x[i+1]-x[i];
			coef[k][0+i*4] = (s2[i+1]-s2[i])/(6*ha);
			coef[k][1+i*4] = s2[i]/2;
			coef[k][2+i*4] = (y[k*n+i+1]-y[k*n+i])/ha-(s2[i+1]+2*s2[i])*(ha/6);
			coef[k][3+i*4] = y[k*n+i];
		}
	}
}

void calcInterfacesZcoord(	float *zi, /* Interfaces depth coordinates */
				int nint, /* Number of interfaces */
				float xs, /* x coordinate */
				int si, /* Spline index */
				float **coef /* Cubic spline coefficients */)
/*< Calculate depth coordinates of the interfaces
 * Note: This function calculates interfaces depth coordinates and stores it
 * in the zi vector.
  >*/
{
	int i; // Loop counter

	for(i=0;i<nint;i++){
		zi[i] = coef[i][si*4+0]*xs*xs*xs+
			coef[i][si*4+1]*xs*xs+
			coef[i][si*4+2]*xs+
			coef[i][si*4+3];
	}
}

float calculateLocationMisfit( float **s, /* NIP sources location */
			   	float *sz, /* Depth coordinates of interfaces */
			   	int nsz, /* NIP sources number for each interface */
			   	float osz, /* sz origin */
			   	float dsz, /* sz sampling */
				int nshot, /* Dimension of the sz vector */
				int itf /* Interface to invert */)
/*< Calculate misfit between NIP sources and interfaces
Note: This function calculates L2 norm of the distances between NIP
sources location in Z and interfaces. The assumption is that NIP sources
are located in the interfaces, so the best velocity model minimize the distance
between then

 >*/
{

	int i; // loop counter
	int l=0; // Splines index
	float *zi; // Temporary vector to store depth coordinates
	float *x; // X coordinates of interface being inverted
	float** coef; // Cubic splines coefficients
	float misfit = 0.; // Misfit sum
	float* szz; // Z coordinates of interface being inverted

	x = sf_floatalloc(nsz);
	szz = sf_floatalloc(nsz);

	for(i=0;i<nsz;i++){
		x[i] = i*dsz+osz;
		szz[i]=sz[i+(itf*nsz)];
	}

	/* Calculate coefficients matrix (interfaces interpolation) */
	coef = sf_floatalloc2(4*(nsz-1),1);
	calculateSplineCoeficients(nsz,x,szz,coef,1);

	zi = sf_floatalloc(1);

	/* Calculate interfaces z coordinates and misfit */
	for(i=0;i<nshot;i++){

		l = (int) (s[i][1]-osz)/dsz;

		calcInterfacesZcoord(zi,1,s[i][1]-x[l],l,coef);

		misfit += fabs(zi[0]-s[i][0]);
	}

	return sqrt(misfit*misfit);
}

float splinevellayer(float xs, /* x coordinate */
		     int si, /* spline index */
		     float** coef /* Coefficients matrix */)
/*< Second layer velocity interpolation using cubic splines
Note: This function returns the velocity value for x coordinate
>*/
{
	return (coef[0][si*4+0]*xs*xs*xs+
		coef[0][si*4+1]*xs*xs+
		coef[0][si*4+2]*xs+
		coef[0][si*4+3]);
}

void updateVelocityModelLateralVariation(
		int *n, /* Velocity model dimension n1=n[0] n2=n[1] */
		float *o, /* Velocity model axis origin o1=o[0] o2=o[1] */
		float *d, /* Velocity model sampling d1=d[0] d2=d[1] */
		float *sv, /* Velocity model disturbance */
		int nsv, /* Dimension of sv the vector */
		float *sz, /* Depth coordinates of interfaces */
		int nsz, /* Dimension sz the vector */
		float osz, /* sz vector origin */
		float dsz, /* sz vector sampling */
		float *vel, /* Velocity model */
		int nvel, /* Dimension of the vel vector */
		float *svx)
/*< Velocity model update, Lateral velocity variation in second layer
Note: This function uses a sv (layers velocity) vector and sz (depth interfaces
coordinates) vector to build the depth velocity model. There is nsv constant
velocity layers (except the second one) in the model and nsv-1 interfaces
separating them.
These interfaces are described with nsz control points in the sz vector and
they are interpolated using natural cubic spline interpolation.
 >*/
{

	int i, j; // Loop counters
	int k; // Layers index
	int l=0; // Splines index
	float z; // Depth coordinate
	float *zi; // Temporary vector to store depth coordinates
	int nx=nsz/(nsv-1); // Number of nodepoints for each interface
	float *x; // x coordinates of the interfaces nodepoints
	float** coef; // interface spline cubic coefficients
	float** coefsx; // second layer's velocity spline cubic coefficients
	float xx; // x coordinates of the velocity model

	x = sf_floatalloc(nx);
	for(i=0;i<nx;i++)
		x[i] = i*dsz+osz;

	/* Calculate coefficients matrix (interfaces interpolation) */
	coef = sf_floatalloc2(4*(nx-1),nsv-1);
	calculateSplineCoeficients(nx,x,sz,coef,nsv-1);

	/* Calculate coefficients matrix (second layer velocity) */
	coefsx = sf_floatalloc2(4*(nx-1),1);
	calculateSplineCoeficients(nx,x,svx,coefsx,1);

	zi = sf_floatalloc(nsv);
	zi[nsv-1] = (n[0]-1)*d[0]+o[0];

	/* Calculate velocity function */
        for(j=0;j<n[1];j++){

		xx = d[1]*j+o[1];
		if(xx>x[l+1]) l++;
		/* Calculate interfaces z coordinates */
		calcInterfacesZcoord(zi,nsv-1,xx-x[l],l,coef);
		k=0;
                for(i=0;i<n[0];i++){
			z = i*d[0]+o[0];
			if(z>zi[k]) k++; // If second layer use velocity function
			vel[(n[0]*j)+i] = (k!=1)? sv[k]:splinevellayer(xx-x[l],l,coefsx);
			//vel[(n[0]*j)+i] = sv[k];
                } /* Loop over depth */

	} /* Loop over distance */
}

void buildSlownessModelFromVelocityModel(int *n, /* Velocity model dimension n1=n[0] n2=n[1] */
			 		 float *o, /* Velocity model axis origin o1=o[0] o2=o[1] */
					 float *d, /* Velocity model sampling d1=d[0] d2=d[1] */
					 float *sv, /* Velociy disturbance */
					 int nsv, /* Dimension of sv vector */
					 float *sz, /* Depth coordinates of interfaces */
					 int nsz, /* Dimension of sz vector */
					 float osz,
					 float dsz,
					 float *vel, /* Velocity model */
					 int nslow, /* Dimension of vel vector */
					 float* svx)
/*< Slowness model build from velocity model
Note: This function is a function wrapper to updateVelocityModel function.
It calls that function to update the velocity model and build the slowness
model matrix using the slowness definition slow=(1.0/(v*v)). 
 >*/
{

	int i, nm; // Loop counters and indexes

	nm =n[0]*n[1];
	updateVelocityModelLateralVariation(n,o,d,sv,nsv,sz,nsz,osz,dsz,vel,nm,svx);

	/* transform velocity to slowness */
	for(i=0;i<nm;i++){
			vel[i] = 1.0/(vel[i]*vel[i]);
	}
}

