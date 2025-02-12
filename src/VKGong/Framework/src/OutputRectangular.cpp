/*
 * VK-Gong code 
 *
 * Author of this version: Àngels Aragonès
 *
 * Based on VK-Gong for Matlab and NeSS Framework Code (Copyright (c) The University of Edinburgh, 2014-2015. All rights reserved. Author: James Perry)
 *
 */

#include "OutputRectangular.h"
#include "Logger.h"
#include "GlobalSettings.h"



#include <cmath>
using namespace std;

OutputRectangular::OutputRectangular(RectangularImperfectPlate *comp, double pan, double x, double y)
{
    logMessage(1, "Creating modal output from %s: %f,%f", comp->getName().c_str(),
	       x, y);

    component = comp;
    this->pan = pan;

    int NF = GlobalSettings::getInstance()->getNumTimesteps();
    data = new double[NF];

    memset(data, 0, NF * sizeof(double));

    // calculate read point vector
    Nphi = comp->getStateSize();
    h = comp->geth();
    rp = new double[Nphi];
    double *ov = comp->getOmega();
    int *kx = comp->getKx();
    int *ky = comp->getKy();
    double Lx = comp->getLx();
    double Ly = comp->getLy();
    double op1 = x ;
    double op2 = y ;
    int i;
    for (i = 0; i < Nphi; i++) {
	rp[i] = sin((double)kx[i] * M_PI * op1 / Lx) * sin((double)ky[i] * M_PI * op2 / Ly);
    }


    logMessage(1, "Adding modal output from %s at position %f,%f, pan %f", comp->getName().c_str(),
	       x, y, pan);


}

OutputRectangular::~OutputRectangular()
{
    delete[] rp;
}

void OutputRectangular::runTimestep(int n)
{

    double *q = component->getU();
    int i;
    double val = 0.0;
    for (i = 0; i < Nphi; i++) {
	val += rp[i] * q[i];

    }
    //data[n] = val / h;
    data[n] = val;
}


