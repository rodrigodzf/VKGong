/*
 * VK-Gong code 
 *
 * Author of this version: Àngels Aragonès
 *
 * Based on VK-Gong for Matlab and NeSS Framework Code (Copyright (c) The University of Edinburgh, 2014-2015. All rights reserved. Author: James Perry)
 *
 */

#include "PhysicalSystem.h"
#include "Logger.h"
#include "WavWriter.h"
#include "GlobalSettings.h"
#include "SettingsManager.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
using namespace std;

PhysicalSystem::PhysicalSystem()
{
    profiler = NULL;

}

PhysicalSystem::~PhysicalSystem()
{
    int i;
    
    for (i = 0; i < components.size(); i++) {	
	delete components.at(i);	
    }
    
    for (i = 0; i < outputs.size(); i++) {
	delete outputs.at(i);
    }

    if (profiler) {
	delete profiler;
    }

    
}

void PhysicalSystem::runTimestepVerlet(int n)
{
    int i;

    if (!profiler) {
	profiler = new Profiler(components.size());
    }


	// run all components on single thread
	for (i = 0; i < components.size(); i++) {
	    profiler->start(i);
	    components.at(i)->runTimestepVerlet(n);
	    profiler->end(i);
	}

}


void PhysicalSystem::runTimestepECS(int n)
{
    int i;

    if (!profiler) {
	profiler = new Profiler(components.size());
    }

	// run all components on single thread
	for (i = 0; i < components.size(); i++) {
	    profiler->start(i);
	    components.at(i)->runTimestepECS(n);
	    profiler->end(i);
	}

}

void PhysicalSystem::endTimestepVerlet(int n)
{
    int i;
    for (i = 0; i < outputs.size(); i++) {
	outputs.at(i)->runTimestep(n);
    }

    for (i = 0; i < components.size(); i++) {
	components.at(i)->swapBuffers(n);
    }
}

void PhysicalSystem::endTimestepECS(int n)
{
    int i;
    for (i = 0; i < outputs.size(); i++) {
	outputs.at(i)->runTimestep(n);
    }

    for (i = 0; i < components.size(); i++) {
	components.at(i)->swapBuffers(n);
    }

    for (i = 0; i < components.size(); i++) {
	components.at(i)->swapBuffersEta(n);
    }
}

double PhysicalSystem::getMaxOutput()
{
    int i, j;
    double max;
    double curr;

    max = 0.0;
    // get absolute maximum across entire output space
    for (i = 0; i < GlobalSettings::getInstance()->getNumTimesteps(); i++) {
	curr = 0.0;
	for (j = 0; j < outputs.size(); j++) {
	    curr += outputs.at(j)->getData()[i];
	}
	curr = fabs(curr);
	if (curr > max) max = curr;
    }
    logMessage(1, "Max output is %f", max);
    return max;
}

void PhysicalSystem::saveOutputs(string outputname, bool individual, bool raw)
{
    int i, j;
    double max;
    GlobalSettings *gs = GlobalSettings::getInstance();
    
    // write raw outputs if requested
    if (raw) {
	logMessage(1, "Saving raw output data");
	for (i = 0; i < outputs.size(); i++) {
	    ostringstream convert;
	    convert << (i+1);
	    string nstr = convert.str();
	    string compname = outputs.at(i)->getComponent()->getName();
	    //string filename = outputname + "-" + compname + "-" + nstr + ".f64";
	    //if (outputs.size() == 1) filename = outputname + ".f64";
	    string filename = "./OutputResults/" + outputname + "-" + compname + "-" + nstr ;
	    if (outputs.size() == 1) filename = "./OutputResults/" + outputname ;
	    outputs.at(i)->saveRawData(filename);
	}
    }

    // normalise the channels if requested
    if (gs->getNormaliseOuts()) {
	logMessage(1, "Normalising output channels");
	for (i = 0; i < outputs.size(); i++) {
	    outputs.at(i)->normalise();
	}
    }

    // always write stereo mix (if there's more than one channel to mix!)
    if (outputs.size() > 1) {
	logMessage(1, "Writing stereo mix output");
	max = getMaxOutput();
	WavWriter mix("./OutputResults/" + outputname + "-mix.wav");
	if (!mix.writeStereoMix(&outputs, max)) {
	    logMessage(5, "Error writing stereo mix output!");
	}
    }

    // write individual output files if requested
    if (individual) {
	logMessage(1, "Writing individual channel wavs");
	for (i = 0; i < outputs.size(); i++) {
	    ostringstream convert;
	    convert << (i+1);
	    string nstr = convert.str();
	    string compname = outputs.at(i)->getComponent()->getName();
	    string filename = "./OutputResults/" + outputname + "-" + compname + "-" + nstr + ".wav";
	    if (outputs.size() == 1) filename = "./OutputResults/" + outputname + ".wav";    	    
	    WavWriter channel(filename);
	    max = outputs.at(i)->getMaxValue();
	    if (!channel.writeMonoWavFile(outputs.at(i), max)) {
		logMessage(5, "Error writing wav file for channel %d!", i+1);
	    }
	}
    }


}

Component *PhysicalSystem::getComponentByName(string name)
{
    int i;
    for (i = 0; i < components.size(); i++) {
	if (name == components.at(i)->getName()) return components.at(i);
    }
    return NULL;
}

