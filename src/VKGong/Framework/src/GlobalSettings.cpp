/*
 * VK-Gong code 
 *
 * Author of this version: Àngels Aragonès
 *
 * Based on VK-Gong for Matlab and NeSS Framework Code (Copyright (c) The University of Edinburgh, 2014-2015. All rights reserved. Author: James Perry)
 *
 */

#include "GlobalSettings.h"
#include "Logger.h"

GlobalSettings::GlobalSettings()
{
    logMessage(1, "GlobalSettings initialising to defaults");

    // sensible default values for everything
    sampleRate = 44100.0;
    duration = 1.0;
    k = 1.0 / sampleRate;

    maxThreads = 1024;
    
    StormerVerlet = true;

    estimate = false;
    maxOut = 1.0;
    
    Nphi_min = 10;
    Npsi_min = 3;
    
    normaliseOuts = false;
    

}

GlobalSettings::~GlobalSettings()
{
}

GlobalSettings *GlobalSettings::instance = NULL;
