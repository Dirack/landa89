/*
	 vfsacrsnh_lib.c (c)
	 
	 Purpose: 'Mvfsacrsnh.c' library.
	 	 
	 Version 1.0
	 
	 Site: https://dirack.github.io
	 
	 Programmer: Rodolfo A. C. Neves (Dirack) 19/09/2019

	 Email:  rodolfo_profissional@hotmail.com

	 License: GPL-3.0 <https://www.gnu.org/licenses/gpl-3.0.txt>.

*/

/*
TODO: Modify macro definition in search window for each interface.
Large windows can make the result oscilate a lot and do not converge
*/
#define MAX_VEL 1.6
#define MIN_VEL 1.45
#define MAX_VEL2 1.8
#define MIN_VEL2 1.6
#define MAX_Z 1.3
#define MIN_Z 0.9
#define APERTURE MAX_VEL-MIN_VEL
#define hMAX 50
#define mMAX 50
#define ITMAX 3000
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <rsf.h>
/*^*/

#define signal(s) ((s<0)?(-1.):(1.))
/*< Signal function >*/
/*^*/

float getRandomNumberBetween0and1(){
/*< Function to get a random number between 0 and 1 >*/

	return (float)(rand()%1000)/1000;
}

float getVfsaIterationTemperature(int iteration,float dampingFactor,float inicialTemperature){
/*< Temperature function for VFSA algorithm >*/

	return inicialTemperature*expf(-dampingFactor*pow(iteration,0.25));

}

/* TODO: Modify this function for multiple interfaces */
void disturbParameters( float temperature, /* Temperature of this interation in VFSA */
			float* disturbedVel, /* Parameters disturbed vector */
			float* originalVel, /* original parameters vector */
			int nv, /*Number of parameters */
			float* disturbedZ, /* Parameters disturbed vector */
			float* originalZ, /* original parameters vector */
			int nz, /*Number of parameters */
			float scale /* Scale to multiply by disturbance */,
			int itf)
/*< Disturb parameters from the previous iteration of VFSA
 Note: It receives a parameter vector and distubs it accordingly to 
VFSA disturb parameters step.
 >*/
{

	float u;
	float disturbance;
	int i;
	int nx=nz/(nv-1);
	// TODO pass max and min values through cmd
	float minz[2]={0.9,1.75};
	float maxz[2]={1.45,1.9};
	float minvel[2]={1.45,1.65};
	float maxvel[2]={1.60,1.85};

	for(i=0;i<nv;i++)		
		disturbedVel[i]=originalVel[i];

	u=getRandomNumberBetween0and1();
				
	disturbance = signal(u - 0.5) * temperature * (pow( (1+temperature),fabs(2*u-1) )-1);

	disturbedVel[itf] = originalVel[itf] + (disturbance*scale*10) * (0.05);

	if (disturbedVel[itf] >= maxvel[itf]) {

		disturbedVel[itf] = maxvel[itf] - (maxvel[itf]-minvel[itf]) * getRandomNumberBetween0and1();
			
	}

	if (disturbedVel[itf] <= minvel[itf]) {

		disturbedVel[itf] = (maxvel[itf]-minvel[itf]) * getRandomNumberBetween0and1() + minvel[itf];
			
	}

	for(i=0;i<nz;i++)
		disturbedZ[i]=originalZ[i];
	
	for(i=(itf*nx);i<(itf*nx+nx);i++){

		u=getRandomNumberBetween0and1();
				
		disturbance = signal(u - 0.5) * temperature * (pow( (1+temperature),fabs(2*u-1) )-1);

		disturbedZ[i] = originalZ[i] + (disturbance*scale);

		if (disturbedZ[i] >= maxz[itf]) {

			disturbedZ[i] = maxz[itf] - (maxz[i]-minz[i]) * getRandomNumberBetween0and1();
			
		}

		if (disturbedZ[i] <= minz[itf]) {

			disturbedZ[i] = (maxz[i]-minz[i]) * getRandomNumberBetween0and1() + minz[itf];
			
		}
	}
}

void nonHyperbolicCRSapp(float t[2*mMAX+1][hMAX], float m0, float dm, float om, float dh, float oh, float t0, float v0, float RN, float RNIP, float BETA){
/*< Non hyperbolic CRS approximation (FOMEL; KAZINNIK, 2013) >*/
	float m0_index=(int)(m0/dm);
	float a1, a2, b2, c1, Fd, Fd1, Fd2;
	int im, ih;
	float m, h, mmh, mph;
	float sinB=sin(BETA),cosB=cos(BETA);
	
	om = om+(m0_index-mMAX)*dm;

	a1=(2*sinB)/(v0);		
	a2=(2*cosB*cosB*t0)/(v0*RN);
	b2=(2*cosB*cosB*t0)/(v0*RNIP);
	c1=2*b2+a1*a1-a2;

	for (im=0; im < 2*mMAX+1; im++){
			
		m=(im*dm+om)-m0;

		for(ih=0;ih<hMAX;ih++){
			
			h=ih*dh+oh;
			mmh=m-h;
			mph=m+h;

			Fd=(t0+a1*m)*(t0+a1*m)+a2*m*m;				
			Fd2=(t0+a1*(mmh))*(t0+a1*(mmh))+a2*(mmh)*(mmh);
			Fd1=(t0+a1*(mph))*(t0+a1*(mph))+a2*(mph)*(mph);					
			t[im][ih]=sqrt((Fd+c1*h*h+sqrt(Fd2*Fd1))*0.5); 
		}
	}

}

float semblance(float m0, float dm, float om, float oh, float dh, float dt, int nt,float t0, float v0,float RN, float RNIP, float BETA, float*** t){
/*< Semblance: Non Hyperbolic CRS approximation with data >*/

	int im, ih, numSamples=0;
	float amplitude=0.;
	float amplitudeSampleSum=0.;
	float amplitudeSquaredSampleSum=0.;
	float semblance=0;
	int tetai;
	float teta[2*mMAX+1][hMAX];
	int m0_index_init, m0_index_end;

	m0_index_init = (int)(m0/dm)-mMAX;
	m0_index_end = (int)(m0/dm)+mMAX;

	nonHyperbolicCRSapp(teta,m0,dm,om,dh,oh,t0,v0,RN,RNIP,BETA);

	for (im=m0_index_init; im < m0_index_end; im++){
			
		for(ih=0;ih<hMAX;ih++){

			tetai=teta[im-m0_index_init][ih]/dt;

			if(tetai>=0 && tetai < nt){
				amplitude = t[im][ih][tetai];
				
				amplitudeSampleSum += amplitude;
					
				amplitudeSquaredSampleSum += (amplitude*amplitude);
			}
				
			numSamples++;
		}
		
	}		

	if(amplitudeSampleSum==0)		
		return semblance=0;
	else
		return semblance=(amplitudeSampleSum*amplitudeSampleSum)/(numSamples*amplitudeSquaredSampleSum);

}

