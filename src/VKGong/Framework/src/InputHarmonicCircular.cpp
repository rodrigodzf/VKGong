/*
 * VK-Gong code 
 *
 * Author of this version: Àngels Aragonès
 *
 * Based on VK-Gong for Matlab and NeSS Framework Code (Copyright (c) The University of Edinburgh, 2014-2015. All rights reserved. Author: James Perry)
 *
 */

#include "InputHarmonicCircular.h"
#include "Logger.h"
#include "GlobalSettings.h"

InputHarmonicCircular::InputHarmonicCircular(CircularImperfectPlate *comp, double th, double r,
				   double startTime, double Amplitude, double freq, double phase, double duration)
    : Input(comp)
{
    logMessage(1, "Creating  harmonic: %f, %f, %f, %f, %f", th, r, startTime,
	       duration, Amplitude);

    
    double tnd = comp->gettnd();
    double h = comp->geth();
    double Rd = comp->getRd();
    double Young = comp->getYoung();
    double e = comp->gete();
    char BC = comp->getBoundaryConditions();
    
    // convert times to discrete timesteps
    this->startTime = (int)floor(timeToTimestep(startTime));   
    this->duration = (int)floor(timeToTimestep(duration));
    
    this->period = GlobalSettings::getInstance()->getSampleRate() / freq;
	this->phase = phase;
    

    if (this->startTime < firstInputTimestep) firstInputTimestep = this->startTime;

    double amp_coef = Rd*Rd*Rd*Rd / Young / h / h / h / h;
    
    this->Amplitude = Amplitude * amp_coef;
    

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

InputHarmonicCircular::~InputHarmonicCircular()
{
    delete[] P;

}

void InputHarmonicCircular::runTimestep(int n, double *s, double *s1, double *s2)
{

    n -= startTime;
    if ((n >= 0) && (n < duration)) {   
	   
    	
   double val = Amplitude * sin(2*M_PI*n / period + phase); 	
    	
    	
	// add it to component's modes
	int i;
	for (i = 0; i < Nphi; i++) {
	    s[i] += val * P[i];

	}

    }
}



