/*
 * VK-Gong code 
 *
 * Author of this version: Àngels Aragonès
 *
 * Based on VK-Gong for Matlab and NeSS Framework Code (Copyright (c) The University of Edinburgh, 2014-2015. All rights reserved. Author: James Perry)
 *
 * Singleton class to hold global settings used by entire simulation
 */
#ifndef _GLOBALSETTINGS_H_
#define _GLOBALSETTINGS_H_

#include <cstdio>
using namespace std;

class GlobalSettings {
 public:
    // get singleton instance
    static GlobalSettings *getInstance() { 
	if (instance == NULL) instance = new GlobalSettings();
	return instance;
    }

    void setSampleRate(double sr) { 
	sampleRate = sr;
	// k is inverse of sampleRate
	k = 1.0 / sampleRate;
    }
    void setDuration(double d) { duration = d; }
    void setK(double kk) { k = kk; }

    void setMaxThreads(int mt) { maxThreads = mt; }

    void setEstimate(bool e) { estimate = e; }

    void setMaxOut(double mo) { maxOut = mo; }
    
    void setStormerVerletScheme(bool sv) { StormerVerlet = sv; }


    void setNormaliseOuts(bool no) { normaliseOuts = no; }


    double getSampleRate() { return sampleRate; }
    double getDuration() { return duration; }
    double getK() { return k; }

    int getNumTimesteps() {
	return (int)(sampleRate * duration);
    }



    int getMaxThreads() { return maxThreads; }

    bool getEstimate() { return estimate; }

    double getMaxOut() { return maxOut; }

    bool getStormerVerletScheme() {return StormerVerlet;}
    
    bool getNormaliseOuts() { return normaliseOuts; }
    
    int getNphiMin() { return Nphi_min;}
    int getNpsiMin() { return Npsi_min;}
    
    int setNphiMin(int Npm) { Nphi_min = Npm;}
    int setNpsiMin(int Npm) { Npsi_min = Npm;}
    

 private:
    static GlobalSettings *instance;
    GlobalSettings();
    ~GlobalSettings();

    double sampleRate;
    double duration;
    double k;

    int maxThreads;

    bool estimate;

    double maxOut;
    
    bool StormerVerlet;

    bool normaliseOuts;
    
    int Nphi_min;
    int Npsi_min;
 
};

#endif
