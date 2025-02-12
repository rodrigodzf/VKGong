/*
 * VK-Gong code 
 *
 * Author of this version: Àngels Aragonès
 *
 * Based on VK-Gong for Matlab and NeSS Framework Code (Copyright (c) The University of Edinburgh, 2014-2015. All rights reserved. Author: James Perry)
 *
 * Represents an output from a rectangular plate component
 */
#ifndef _OUTPUT_RECTANGULAR_H_
#define _OUTPUT_RECTANGULAR_H_

#include "Output.h"
#include "RectangularImperfectPlate.h"

class OutputRectangular : public Output {
 public:
    OutputRectangular(RectangularImperfectPlate *comp, double pan, double x, double y);
    virtual ~OutputRectangular();

    virtual void runTimestep(int n);

 protected:
    double *rp;
    int Nphi;
    double h;
};


#endif
