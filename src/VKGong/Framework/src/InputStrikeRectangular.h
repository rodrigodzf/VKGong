/*
 * VK-Gong code 
 *
 * Author of this version: Àngels Aragonès
 *
 * Based on VK-Gong for Matlab and NeSS Framework Code (Copyright (c) The University of Edinburgh, 2014-2015. All rights reserved. Author: James Perry)
 *
 * Represents a strike on a rectangular plate
 */

#ifndef _INPUT_STRIKE_RECTANGULAR_H_
#define _INPUT_STRIKE_RECTANGULAR_H_

#include "Input.h"
#include "RectangularImperfectPlate.h"

#include <gsl/gsl_sf_bessel.h>


class InputStrikeRectangular : public Input {
 public:
    InputStrikeRectangular(RectangularImperfectPlate *comp, double x, double y, double startTime, double duration, double fm);
    virtual ~InputStrikeRectangular();
    
    virtual void runTimestep(int n, double *s, double *s1, double *s2);

 protected:
    double *P;
    int Nphi;
    double fm;
    

};


#endif
