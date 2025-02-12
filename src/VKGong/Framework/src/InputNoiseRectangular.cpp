/*
 * VK-Gong code 
 *
 * Author of this version: Àngels Aragonès
 *
 * Based on VK-Gong for Matlab and NeSS Framework Code (Copyright (c) The University of Edinburgh, 2014-2015. All rights reserved. Author: James Perry)
 *
 */

#include "InputNoiseRectangular.h"
#include "Logger.h"
#include "GlobalSettings.h"

#include <cstdlib>

InputNoiseRectangular::InputNoiseRectangular(RectangularImperfectPlate *comp, double x, double y,
				   double startTime, double Amplitude, double fmin, double fmax, double deltaf, double duration)
    : Input(comp)
{
    logMessage(1, "Creating  noise: %f, %f, %f, %f, %f", x, y, startTime,
	       duration, Amplitude);

    char BC = comp->getBoundaryConditions();
    
    // convert times to discrete timesteps
    this->startTime = (int)floor(timeToTimestep(startTime));   
    this->duration = (int)floor(timeToTimestep(duration));

    if (this->startTime < firstInputTimestep) firstInputTimestep = this->startTime;    


    // generate P
   Nphi = comp->getStateSize();
   P = new double[Nphi];
   double Lx = comp->getLx();
   double Ly = comp->getLy();
   double h = comp->geth();
   double rho = comp->getRho();
   double *Cnorm = comp->getCnorm();
   double *ov = comp->getOmega();
   int *kx = comp->getKx();
   int *ky = comp->getKy();

   
   int i;       

    if (BC == 's') {
		for (i = 0; i < Nphi; i++) {	
			P[i] = sin(x * M_PI * kx[i]/Lx) * sin(y * M_PI * ky[i]/Ly);
			P[i] = ((P[i]/rho)/(Lx*Ly/4.0))/h;
			P[i] /= Cnorm[i];

		}
    }
    

}

InputNoiseRectangular::~InputNoiseRectangular()
{
    delete[] P;

}

void InputNoiseRectangular::runTimestep(int n, double *s, double *s1, double *s2)
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



