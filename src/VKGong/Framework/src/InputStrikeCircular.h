/*
 * VK-Gong code 
 *
 * Author of this version: Àngels Aragonès
 *
 * Based on VK-Gong for Matlab and NeSS Framework Code (Copyright (c) The University of Edinburgh, 2014-2015. All rights reserved. Author: James Perry)
 *
 * Represents a strike on a circular plate
 */

#ifndef _INPUT_STRIKE_CIRCULAR_H_
#define _INPUT_STRIKE_CIRCULAR_H_

#include "Input.h"
#include "CircularImperfectPlate.h"

#include <gsl/gsl_sf_bessel.h>


class InputStrikeCircular : public Input {
 public:
    InputStrikeCircular(CircularImperfectPlate *comp, double th, double r, double startTime, double duration, double fm);
    virtual ~InputStrikeCircular();
    
    virtual void runTimestep(int n, double *s, double *s1, double *s2);


 protected:
    double *P;
    int Nphi;
    double fm;
    
    

};


#endif
