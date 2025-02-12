/*
 * VK-Gong code 
 *
 * Author of this version: Àngels Aragonès
 *
 * Based on VK-Gong for Matlab and NeSS Framework Code (Copyright (c) The University of Edinburgh, 2014-2015. All rights reserved. Author: James Perry)
 *
 * Represents an output from a circular plate component
 */
#ifndef _OUTPUT_CIRCULAR_H_
#define _OUTPUT_CIRCULAR_H_

#include "Output.h"
#include "CircularImperfectPlate.h"

class OutputCircular : public Output {
 public:
    OutputCircular(CircularImperfectPlate *comp, double pan, double r, double theta);
    virtual ~OutputCircular();

    virtual void runTimestep(int n);
  
 protected:
    double *rp;
    int Nphi;
    double hd, nu;
  
};

#endif
