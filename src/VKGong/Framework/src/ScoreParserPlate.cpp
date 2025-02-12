/*
 * VK-Gong code 
 *
 * Author of this version: Àngels Aragonès
 *
 * Based on VK-Gong for Matlab and NeSS Framework Code (Copyright (c) The University of Edinburgh, 2014-2015. All rights reserved. Author: James Perry)
 *
 */

#include "ScoreParserPlate.h"
#include "Logger.h"
#include "GlobalSettings.h"
#include "CircularImperfectPlate.h"
#include "InputStrikeCircular.h"
#include "InputHarmonicCircular.h"
#include "InputNoiseCircular.h"
#include "RectangularImperfectPlate.h"
#include "InputStrikeRectangular.h"
#include "InputHarmonicRectangular.h"
#include "InputNoiseRectangular.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cstring>
using namespace std;

ScoreParserPlate::ScoreParserPlate(string filename) : ScoreParser(filename)
{
}

ScoreParserPlate::~ScoreParserPlate()
{
}

bool ScoreParserPlate::parse(PhysicalSystem *physicalSystem)
{
    durationOnly = false;
    this->physicalSystem = physicalSystem;
    return parseTextFile();
}

Input *ScoreParserPlate::parseStrike(istream &in)
{
    double startTime, thx, ry, duration, fm;
    string compname;
    Input *strike;

    in >> startTime >> compname >> thx >> ry >> duration >> fm;
        
    if (in.fail()) {
	logMessage(3, "ScoreParserPlate: error parsing strike definition");
	return NULL;
    }

    Component *comp = physicalSystem->getComponentByName(compname);
    if (comp == NULL) {
	logMessage(3, "ScoreParserPlate: unrecognised component name '%s' in strike",
		   compname.c_str());
	return NULL;
    }

	CircularImperfectPlate *cIp = dynamic_cast<CircularImperfectPlate*>(comp);
	if (cIp == NULL){
		RectangularImperfectPlate *rIp = dynamic_cast<RectangularImperfectPlate*>(comp);
		if (rIp == NULL){
			logMessage(3, "ScoreParserPlate: error parsing component");
		}
		else{
			strike = new InputStrikeRectangular(rIp, thx, ry, startTime, duration, fm); 
		}
	}
	else{

		strike = new InputStrikeCircular(cIp, thx, ry, startTime, duration, fm); 
	}

    
    comp->addInput(strike);

    
    return strike;
}

Input *ScoreParserPlate::parseHarmonicInput(istream &in)
{
	double startTime, thx, ry, Amplitude, freq, phase, duration;
	string compname;
	Input *harmonic;

	in >> startTime >> compname >> thx >> ry >> Amplitude >> freq >> phase >> duration;
	if (in.fail()) {
	logMessage(3, "ScoreParserPlate: error parsing harmonic input definition");
	return NULL;
	}

	Component *comp = physicalSystem->getComponentByName(compname);
	if (comp == NULL) {
	logMessage(3, "ScoreParserPlate: unrecognised component name '%s' in harmonic input",
		   compname.c_str());
	return NULL;
	}




	CircularImperfectPlate *cIp = dynamic_cast<CircularImperfectPlate*>(comp);
	if (cIp == NULL){
		RectangularImperfectPlate *rIp = dynamic_cast<RectangularImperfectPlate*>(comp);
		if (rIp == NULL){
			logMessage(3, "ScoreParserPlate: error parsing component");
		}
		else{
			harmonic = new InputHarmonicRectangular(rIp, thx, ry, startTime, Amplitude, freq, phase, duration); 
		}
	}
	else{
		harmonic = new InputHarmonicCircular(cIp, thx, ry, startTime, Amplitude, freq, phase, duration); 
	}

	
	comp->addInput(harmonic);
	return harmonic;
	
}

Input *ScoreParserPlate::parseNoiseInput(istream &in)
{
	double startTime, thx, ry, Amplitude, fmin, fmax, deltaf, duration;
	string compname;
	Input *noise;

	in >> startTime >> compname >> thx >> ry >> Amplitude >> fmin >> fmax >> deltaf >> duration;
	if (in.fail()) {
	logMessage(3, "ScoreParserPlate: error parsing harmonic input definition");
	return NULL;
	}

	Component *comp = physicalSystem->getComponentByName(compname);
	if (comp == NULL) {
	logMessage(3, "ScoreParserPlate: unrecognised component name '%s' in harmonic input",
		   compname.c_str());
	return NULL;
	}

	CircularImperfectPlate *cIp = dynamic_cast<CircularImperfectPlate*>(comp);
	if (cIp == NULL){
		RectangularImperfectPlate *rIp = dynamic_cast<RectangularImperfectPlate*>(comp);
		if (rIp == NULL){
			logMessage(3, "ScoreParserPlate: error parsing component");
		}
		else{
			noise = new InputNoiseRectangular(rIp, thx, ry, startTime, Amplitude, fmin, fmax, deltaf, duration); 
		}
	}
	else{
		noise = new InputNoiseCircular(cIp, thx, ry, startTime, Amplitude, fmin, fmax, deltaf, duration); 
	}

	
	comp->addInput(noise);
	return noise;
}

double ScoreParserPlate::getDuration()
{
    durationOnly = true;
    duration = -1.0;
    parseTextFile();
    return duration;
}

int ScoreParserPlate::handleItem(string type, istream &in)
{
    int result = 1;

    if (durationOnly) {
	if (type == "duration") {
	    in >> duration;
	    logMessage(3, "Setting duration to %f", duration);
	    return 1;
	}
    }
    else {
	if (type == "strike") {
	    Input *strike = parseStrike(in);
	    if (!strike) {	    	
		result = 0;
	    }

	}
	else if (type == "harmonic") {
	    Input *harmonic = parseHarmonicInput(in);
	    if (!harmonic) {
		result = 0;
	    }
	}
	else if (type == "noise") {
	    Input *noise = parseNoiseInput(in);
	    if (!noise) {
		result = 0;
	    }
	}
	
	else if (type == "duration") {
	    // already done
	}
	else {
	    logMessage(3, "ScoreParserPlate: unrecognised line '%s' in score file",
		       type.c_str());
	    result = 0;
	}
    }


    return result;
}
