/*
 * VK-Gong code 
 *
 * Author of this version: Àngels Aragonès
 *
 * Based on VK-Gong for Matlab and NeSS Framework Code (Copyright (c) The University of Edinburgh, 2014-2015. All rights reserved. Author: James Perry)
 *
 * Represents a noise input on a circular plate
 */

#ifndef _INPUT_NOISE_CIRCULAR_H_
#define _INPUT_NOISE_CIRCULAR_H_

#include "Input.h"
#include "CircularImperfectPlate.h"

#include <gsl/gsl_sf_bessel.h>



class InputNoiseCircular : public Input {
 public:
	InputNoiseCircular(CircularImperfectPlate *comp, double th, double r,
					   double startTime, double Amplitude, double fmin, double fmax, double deltaf, double duration);
    virtual ~InputNoiseCircular();
    
    virtual void runTimestep(int n, double *s, double *s1, double *s2);


 protected:
    double *P;
    int Nphi;
    double Amplitude;
    double *RandPhase;
    double Tmin;
    double Tmax;
    double deltaT;
    int Nt;
    
    

};


#endif
