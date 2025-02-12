/*
 * VK-Gong code 
 *
 * Author of this version: Àngels Aragonès
 *
 * Based on VK-Gong for Matlab and NeSS Framework Code (Copyright (c) The University of Edinburgh, 2014-2015. All rights reserved. Author: James Perry)
 *
 */

#include "InputHarmonicRectangular.h"
#include "Logger.h"
#include "GlobalSettings.h"

InputHarmonicRectangular::InputHarmonicRectangular(RectangularImperfectPlate *comp, double x, double y,
				   double startTime, double Amplitude, double freq, double phase, double duration)
    : Input(comp)
{
    logMessage(1, "Creating  harmonic: %f, %f, %f, %f, %f", x, y, startTime,
	       duration, Amplitude);

    
    
    char BC = comp->getBoundaryConditions();
    
    // convert times to discrete timesteps
    this->startTime = (int)floor(timeToTimestep(startTime));   
    this->duration = (int)floor(timeToTimestep(duration));
    
    this->period = GlobalSettings::getInstance()->getSampleRate() / freq;
	this->phase = phase;
    
    if (this->startTime < firstInputTimestep) firstInputTimestep = this->startTime;    
    
    this->Amplitude = Amplitude;
    

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

    if (BC == 's') {
		for ( int i = 0; i < Nphi; i++) {	
			P[i] = sin(x * M_PI * kx[i]/Lx) * sin(y * M_PI * ky[i]/Ly);
			P[i] = ((P[i]/rho)/(Lx*Ly/4.0))/h;
			P[i] /= Cnorm[i];
		}
	}
}

InputHarmonicRectangular::~InputHarmonicRectangular()
{
    delete[] P;

}

void InputHarmonicRectangular::runTimestep(int n, double *s, double *s1, double *s2)
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



