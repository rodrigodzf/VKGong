/*
 * VK-Gong code 
 *
 * Author of this version: Àngels Aragonès
 *
 * Based on VK-Gong for Matlab and NeSS Framework Code (Copyright (c) The University of Edinburgh, 2014-2015. All rights reserved. Author: James Perry)
 *
 */

#include "PhysicalSystemParserCircular.h"
#include "Logger.h"
#include "GlobalSettings.h"

#include <iostream>
#include <fstream>
#include <sstream>
using namespace std;

PhysicalSystemParserCircular::PhysicalSystemParserCircular(string filename)
    : PhysicalSystemParser(filename)
{
}

PhysicalSystemParserCircular::~PhysicalSystemParserCircular()
{
}

PhysicalSystem *PhysicalSystemParserCircular::parse()
{
    logMessage(1, "PhysicalSystemParserCircular: attempting to parse %s", filename.c_str());

    physicalSystem = new PhysicalSystem();

    if (!parseTextFile()) {       
	delete physicalSystem;
	physicalSystem = NULL;
    }
    
    return physicalSystem;
}

int PhysicalSystemParserCircular::handleItem(string type, istream &in)
{
    int result = 1;
    double fs = 0.0;
    int Nphimin, Npsimin;
    
    GlobalSettings *gs = GlobalSettings::getInstance();
    

    if (type == "CircularImperfectPlate") {    
    	

   
	CircularImperfectPlate *mp = readCircularImperfectPlate(in);    

 
	if (mp) {	
	    physicalSystem->addComponent(mp);
	}
	else {
	    result = 0;
	}
    }
    else if (type == "output") {    	

    	OutputCircular *op = readOutput(in);

     
	if (op) {
	    physicalSystem->addOutput(op);
	}
	else {
	    result = 0;
	}
    }
    else if (type == "samplerate") {
		in >> fs;
		
		gs->setSampleRate(fs);
		
		printf("Sample rate introduced by user: %f \n", fs);
	
    }
    else if (type == "HfileDim"){
    	
    	in >> Nphimin >> Npsimin;
    	
    	gs->setNphiMin(Nphimin);
    	gs->setNpsiMin(Npsimin);
    	
    	printf("H file dimensions introduced by user: Nphi = %d Npsi = %d \n", Nphimin, Npsimin);
    	
    }
    else {
	logMessage(3, "PhysicalSystemParserCircular: unrecognised line '%s'", type.c_str());
	result = 0;
    }
    return result;
}

CircularImperfectPlate *PhysicalSystemParserCircular::readCircularImperfectPlate(istream &in)
{
    string name;
    int Nphi, Npsi, modeType, Nr, Nth, fs;
    double Rd, h, H, nu, Young, rho, dFac, dExp, dCons, tau2, KR, KT;
    char BC, ImperfectionType;

    GlobalSettings *gs = GlobalSettings::getInstance();    
    
    in >> name >> nu >> Young >> rho >> Rd >> h >> H >> ImperfectionType >> tau2 >> modeType >> dFac >> dExp >> dCons >> BC >> KR >> KT >> Nphi >> Npsi >> Nr >> Nth;
    
    fs = (int)gs->getSampleRate();
    
    if (in.fail()) {
	logMessage(3, "PhysicalSystemParserCircular: error parsing plate definition in physicalSystem file");
	return NULL;
    }    
    
    
    logMessage(3, "Creating circular plate %s, %d, %d, %f, %f, %f, %c, %f, %f, %f, %f, %c, %d, %d", 
	       name.c_str(), Nphi, Npsi, Rd, h, H, ImperfectionType, modeType, nu, Young, rho, BC, Nr, Nth);

    
    CircularImperfectPlate *mp = new CircularImperfectPlate(name, Nphi, Npsi, Rd, h, H, ImperfectionType, tau2, modeType,
 	       nu, Young, rho, BC, KR, KT, Nr, Nth, dFac, dExp, dCons, fs);
    
    
    return mp;
}

OutputCircular *PhysicalSystemParserCircular::readOutput(istream &in)
{
    string compname;
    double th, r, pan = 0.0;

    in >> compname >> th >> r;
    if (in.fail()) {
	logMessage(3, "PhysicalSystemParserCircular: error parsing output definition");
	return NULL;
    }
    Component *comp = physicalSystem->getComponentByName(compname);

    if (comp == NULL) {
	logMessage(3, "PhysicalSystemParserCircular: unrecognised component name '%s'",
		   compname.c_str());
	return NULL;
    }
	    
    CircularImperfectPlate *mp = dynamic_cast<CircularImperfectPlate*>(comp);

    OutputCircular *op = new OutputCircular(mp, pan, r, th);    

    logMessage(3, "Creating output %s, %f, %f, %f", compname.c_str(), th, r, pan);
    return op;
}
