/*
 * VK-Gong code 
 *
 * Author of this version: Àngels Aragonès
 *
 * Based on VK-Gong for Matlab and NeSS Framework Code (Copyright (c) The University of Edinburgh, 2014-2015. All rights reserved. Author: James Perry)
 *
 */

#include "OutputCircular.h"
#include "Logger.h"
#include "GlobalSettings.h"

#include <gsl/gsl_sf_bessel.h>
#include <cmath>

using namespace std;

OutputCircular::OutputCircular(CircularImperfectPlate *comp, double pan, double r, double theta)
{
    logMessage(1, "Creating modal output from %s: %f,%f", comp->getName().c_str(),
	       r, theta);

    component = comp;
    this->pan = pan;

    int NF = GlobalSettings::getInstance()->getNumTimesteps();
    data = new double[NF];

    memset(data, 0, NF * sizeof(double));

    // calculate read point vector
    Nphi = comp->getStateSize();
    hd = comp->geth();
    nu = comp->getnu();
    rp = new double[Nphi];
    double *ov = comp->getOmega();
    double Kkn;
    int *k = comp->getk();
    int *c = comp->getc();
    double *xkn = comp->getXi();
    int Nr = comp->getNr();
    char BC = comp->getBoundaryConditions();
    double KR = comp->getKR();

    double J0,J1,J2,I0,I1,I2,J,I, JJ0, II0;
    double Jtild,Itild, Jkn, Ikn;

    int i;

    
    
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
				rp[i] = (JJ0 - (Jtild*II0/Itild))*cos(k[i]*theta);
	
			}else{
				rp[i] = (JJ0 - (Jtild*II0/Itild))*sin(k[i]*theta);
	
			}

			// Normalisation
			
			Kkn  = comp->Norm_Modes(k[i],xkn[i], 1);
			rp[i] = rp[i] * Kkn;      
			
		}
    }
    else{
    	for (i = 0; i < Nphi; i++) {
    		
    		Jkn = gsl_sf_bessel_Jn(k[i],xkn[i]);
			Ikn = gsl_sf_bessel_In(k[i],xkn[i]);
			
			J = gsl_sf_bessel_Jn(k[i],xkn[i]*r);
			I = gsl_sf_bessel_In(k[i],xkn[i]*r);
			
			if (c[i] == 1){
				rp[i] = (Ikn*J - Jkn*I)*cos(k[i]*theta);
	
			}else{
				rp[i] = (Ikn*J - Jkn*I)*sin(k[i]*theta);
	
			}
	
			// Normalisation
	
			
			Kkn  = comp->Norm_Modes(k[i],xkn[i], 1);
			
			rp[i] = rp[i] * Kkn;        
			
		}
    }
        
    logMessage(1, "Adding modal output from %s at position %f,%f, pan %f", comp->getName().c_str(),
	       theta, r, pan);


}



OutputCircular::~OutputCircular()
{
    delete[] rp;
}

void OutputCircular::runTimestep(int n)
{

    double *q = component->getU();
    int i;
    double val = 0.0;
    for (i = 0; i < Nphi; i++) {
	val += rp[i] * q[i];	
	}
   
    data[n] = val; 
    


}




