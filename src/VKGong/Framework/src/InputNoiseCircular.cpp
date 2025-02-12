/*
 * VK-Gong code 
 *
 * Author of this version: Àngels Aragonès
 *
 * Based on VK-Gong for Matlab and NeSS Framework Code (Copyright (c) The University of Edinburgh, 2014-2015. All rights reserved. Author: James Perry)
 *
 */

#include "InputNoiseCircular.h"
#include "Logger.h"
#include "GlobalSettings.h"

#include <cstdlib>

InputNoiseCircular::InputNoiseCircular(CircularImperfectPlate *comp, double th, double r,
				   double startTime, double Amplitude, double fmin, double fmax, double deltaf, double duration)
    : Input(comp)
{
    logMessage(1, "Creating  noise: %f, %f, %f, %f, %f", th, r, startTime,
	       duration, Amplitude);

    // Non-dimensionalization is included in the noise
    double tnd = comp->gettnd();
    double h = comp->geth();
    double Rd = comp->getRd();
    double Young = comp->getYoung();
    double e = comp->gete();
    char BC = comp->getBoundaryConditions();
    
    // convert times to discrete timesteps
    this->startTime = (int)floor(timeToTimestep(startTime));   
    this->duration = (int)floor(timeToTimestep(duration));

    double fs = GlobalSettings::getInstance()->getSampleRate();

    if (this->startTime < firstInputTimestep) firstInputTimestep = this->startTime;

    double amp_coef = Rd*Rd*Rd*Rd / Young / h / h / h / h;
    
    this->Amplitude= Amplitude * amp_coef;
    
    this->Tmin = fs / fmax;
    this->Tmax = fs / fmin;
    this->deltaT = fs / deltaf;
    this->Nt = (int)((Tmax-Tmin)/deltaT);
    
    // Create a vector of uniformly distributed random phases
    double max = 0;
    
    RandPhase = new double [Nt];
    for (int i=0; i<Nt; i++){
    	RandPhase[i] = (double)rand();
    	if (RandPhase[i] > max) { max = RandPhase[i]; } // This should be improved using <random>
    }
    
    for (int i=0; i<Nt; i++){ RandPhase[i] = RandPhase[i] / max; }
    

    // generate P
    Nphi = comp->getStateSize();
    P = new double[Nphi];
    double *Cnorm = comp->getCnorm();
    double *ov = comp->getOmega();
    double *xkn = comp->getXi();
    int i;

    double J0,J1,J2,I0,I1,I2,J,I, JJ0, II0;
    double Jtild,Itild, Jkn, Ikn;
    double Kkn;
    int *k = comp->getk();
    int *c = comp->getc();
    double KR = comp->getKR();
    
    int Nr = comp->getNr();
    double nu = comp->getnu();
    double dr =0;
    

    dr =  1.0/(Nr-1); // Assuming a normalized radius Rd/Rd = 1;

    if ((BC == 'f') || (BC == 'e')){
		for (i = 0; i < Nphi; i++) {
	
			J2 = gsl_sf_bessel_Jn(k[i]-2,xkn[i]);
			J1 = gsl_sf_bessel_Jn(k[i]-1,xkn[i]);
			J0 = gsl_sf_bessel_Jn(k[i],xkn[i]);
			I2 = gsl_sf_bessel_In(k[i]-2,xkn[i]);
			I1 = gsl_sf_bessel_In(k[i]-1,xkn[i]);
			I0 = gsl_sf_bessel_In(k[i],xkn[i]);
			
			Jtild=ov[i]*J2 + ((nu-2*k[i]+1)*xkn[i] + KR)*J1 + (k[i]*(k[i]+1)*(1-nu) - KR*k[i])*J0;
			Itild=ov[i]*I2 + ((nu-2*k[i]+1)*xkn[i] + KR)*I1 + (k[i]*(k[i]+1)*(1-nu) - KR*k[i])*I0;
	
			JJ0 = gsl_sf_bessel_Jn(k[i],xkn[i]*r);
			II0 = gsl_sf_bessel_In(k[i],xkn[i]*r);
	
			if (c[i] == 1){
				P[i] = (JJ0 - (Jtild*II0/Itild))*cos(k[i]*th);
	
			}else{
				P[i] = (JJ0 - (Jtild*II0/Itild))*sin(k[i]*th);
	
			}
	
			// Normalisation
			Kkn  = comp->Norm_Modes(k[i],xkn[i], 1);
	
			P[i] = P[i] * Kkn * e / Rd / Rd;
			
			P[i] /= Cnorm[i];      

		}
    }
    else {
    	for (i = 0; i < Nphi; i++) {
    		
			Jkn = gsl_sf_bessel_Jn(k[i],xkn[i]);
			Ikn = gsl_sf_bessel_In(k[i],xkn[i]);
			
			J = gsl_sf_bessel_Jn(k[i],xkn[i]*r);
			I = gsl_sf_bessel_In(k[i],xkn[i]*r);
	
	
			if (c[i] == 1){
				P[i] = (Ikn*J - Jkn*I)*cos(k[i]*th);
	
			}else{
				P[i] = (Ikn*J - Jkn*I)*sin(k[i]*th);
	
			}
	
			// Normalisation
			Kkn  = comp->Norm_Modes(k[i],xkn[i], 1);
	
			P[i] = P[i] * Kkn * e / Rd / Rd;
			
			P[i] /= Cnorm[i];      
    		
		}
    }

}

InputNoiseCircular::~InputNoiseCircular()
{
    delete[] P;

}

void InputNoiseCircular::runTimestep(int n, double *s, double *s1, double *s2)
{
	double period =  Tmin;
	double val=0;
	
    n -= startTime;
    if ((n >= 0) && (n < duration)) {  
    	for(int i=0; i < Nt; i++){
    		val = val + cos(2*M_PI*(n / period + RandPhase[i]));
    		period = period + deltaT;
    	}
    	
       	val = Amplitude *val;
    	
	// add it to component's modes
	for (int i = 0; i < Nphi; i++) {
	    s[i] += val * P[i];
	}

    }
}



