/*
 * VK-Gong code 
 *
 * Author of this version: Àngels Aragonès
 *
 * Based on VK-Gong for Matlab and NeSS Framework Code (Copyright (c) The University of Edinburgh, 2014-2015. All rights reserved. Author: James Perry)
 *
 */
#include "Output.h"
#include "GlobalSettings.h"
#include "Logger.h"
#include "SettingsManager.h"

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
using namespace std;

Output::Output()
{
    data = NULL;
}

Output::~Output()
{
    if (data) delete[] data;


}

void Output::runTimestep(int n)
{


}

void Output::saveRawData(string filename)
{
	string filenamebin = filename + ".bin";
	string filenametxt = filename + ".txt";
	string textline;
    ofstream ofbin(filenamebin.c_str(), ios::out | ios::binary);
    if (!ofbin.good()) {
	logMessage(5, "Failed to create raw output file %s", filenamebin.c_str());
	return;
    }
    ofbin.write((const char *)data,
	     GlobalSettings::getInstance()->getNumTimesteps() * sizeof(double));
    ofbin.close();
    ofstream oftxt(filenametxt.c_str());
    if (!oftxt.good()) {
	logMessage(5, "Failed to create raw output file %s", filenametxt.c_str());
	return;
    }
    for (int i = 0; i < GlobalSettings::getInstance()->getNumTimesteps(); i++){
    	oftxt << "data[" << i << "] = " << data[i] << "\n";
    }
    oftxt.close();
    
    
    
    
}

void Output::highPassFilter()
{
    int i;
    // FIXME: the loop is using the value that was just overwritten in previous iteration, is that correct??
    for (i = 1; i < GlobalSettings::getInstance()->getNumTimesteps(); i++) {
	data[i] -= data[i-1];
    }
}

double Output::getMaxValue()
{
    int i;
    double max = 0.0;
    for (i = 0; i < GlobalSettings::getInstance()->getNumTimesteps(); i++) {
	double v = fabs(data[i]);
	if (v > max) max = v;
    }
    return max;
}

void Output::normalise()
{
    double max = getMaxValue();
    int i;
    for (i = 0; i < GlobalSettings::getInstance()->getNumTimesteps(); i++) {
	data[i] /= max;
    }
}




