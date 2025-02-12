/*
 * VK-Gong code 
 *
 * Author of this version: Àngels Aragonès
 *
 * Based on VK-Gong for Matlab and NeSS Framework Code (Copyright (c) The University of Edinburgh, 2014-2015. All rights reserved. Author: James Perry)
 *
 */

#include "InputStrikeRectangular.h"
#include "Logger.h"

InputStrikeRectangular::InputStrikeRectangular(RectangularImperfectPlate *comp, double x, double y,
				   double startTime, double duration,
				   double fm)
    : Input(comp)
{
    logMessage(1, "Creating  strike: %f, %f, %f, %f, %f", x, y, startTime,
	       duration, fm);

   
    char BC = comp->getBoundaryConditions();
    
    // convert times to discrete timesteps
    this->startTime = (int)floor(timeToTimestep(startTime));   
    this->duration = (int)floor(timeToTimestep(duration));

    

    if (this->startTime < firstInputTimestep) firstInputTimestep = this->startTime;

    
    
    this->fm = fm;

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

InputStrikeRectangular::~InputStrikeRectangular()
{
    delete[] P;

}

void InputStrikeRectangular::runTimestep(int n, double *s, double *s1, double *s2)
{

    n -= startTime;
    if ((n >= 0) && (n < (2*duration + 1))) {  // Unity is added to duration to match the results from matlab, since duration is down rounded 
    	// work out value of strike at this point in time
    	
        	
        double val = 0.5 * fm * (1.0 + cos(M_PI * ((double)n - (double)duration) / (double)duration ));
        
        

    	// add it to component's modes
    	int i;
    	for (i = 0; i < Nphi; i++) {
    	    s[i] += val * P[i];   	    
    	   
    	}


    }
}



