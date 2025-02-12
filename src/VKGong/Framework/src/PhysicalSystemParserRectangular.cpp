/*
 * VK-Gong code 
 *
 * Author of this version: Àngels Aragonès
 *
 * Based on VK-Gong for Matlab and NeSS Framework Code (Copyright (c) The University of Edinburgh, 2014-2015. All rights reserved. Author: James Perry)
 *
 */

#include "PhysicalSystemParserRectangular.h"
#include "Logger.h"
#include "GlobalSettings.h"

#include <iostream>
#include <fstream>
#include <sstream>
using namespace std;

PhysicalSystemParserRectangular::PhysicalSystemParserRectangular(string filename)
    : PhysicalSystemParser(filename)
{
}

PhysicalSystemParserRectangular::~PhysicalSystemParserRectangular()
{
}

PhysicalSystem *PhysicalSystemParserRectangular::parse()
{
    logMessage(1, "PhysicalSystemParserRectangular: attempting to parse %s", filename.c_str());

    physicalSystem = new PhysicalSystem();

    if (!parseTextFile()) {       
	delete physicalSystem;
	physicalSystem = NULL;
    }
    
    return physicalSystem;
}

int PhysicalSystemParserRectangular::handleItem(string type, istream &in)
{
    int result = 1;
    double fs = 0.0;
    int Nphimin, Npsimin;
        
    GlobalSettings *gs = GlobalSettings::getInstance();


    if (type == "RectangularImperfectPlate") {    
    	
	RectangularImperfectPlate *mp = readRectangularImperfectPlate(in);    

	if (mp) {	
	    physicalSystem->addComponent(mp);
	}
	else {
	    result = 0;
	}
    }
    else if (type == "output") {    	

    	OutputRectangular *op = readOutput(in);

     
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
	logMessage(3, "PhysicalSystemParserRectangular: unrecognised line '%s'", type.c_str());
	result = 0;
    }
    return result;
}

RectangularImperfectPlate *PhysicalSystemParserRectangular::readRectangularImperfectPlate(istream &in)
{
    string name;
    int Nphi, Npsi, modeType, Nx, Ny, fs;
    double Lx, Ly, h, H, nu, Young, rho, dFac, dExp, dCons, xWidth, yWidth;
    char BC, ImperfectionType;
    
    GlobalSettings *gs = GlobalSettings::getInstance(); 

    in >> name >> nu >> Young >> rho >> Lx >> Ly >> h >> H >> ImperfectionType >> xWidth >> yWidth >> modeType >> dFac >> dExp >> dCons >> Nphi >> Npsi >> BC >> Nx >> Ny;
     
    fs = (int)gs->getSampleRate();
    
    if (in.fail()) {
	logMessage(3, "PhysicalSystemParserRectangular: error parsing plate definition in physicalSystem file");
	return NULL;
    }    
    
    logMessage(3, "Creating rectangular plate %s, %d, %d, %f, %f, %f, %f, %c, %d, %f, %f, %f, %c, %d, %d", 
	       name.c_str(), Nphi, Npsi, Lx, Ly, h, H, ImperfectionType, modeType, nu, Young, rho, BC, Nx, Ny);

    
    RectangularImperfectPlate *mp = new RectangularImperfectPlate(name, Nphi, Npsi, Lx, Ly, h, H, ImperfectionType, xWidth, yWidth, modeType,
 	       nu, Young, rho, BC, Nx, Ny, dFac, dExp, dCons, fs);
    
    
    return mp;
}

OutputRectangular *PhysicalSystemParserRectangular::readOutput(istream &in)
{
    string compname;
    double x, y, pan=0.0;

    in >> compname >> x >> y;
    if (in.fail()) {
	logMessage(3, "PhysicalSystemParserRectangular: error parsing output definition");
	return NULL;
    }
    Component *comp = physicalSystem->getComponentByName(compname);

    if (comp == NULL) {
	logMessage(3, "PhysicalSystemParserRectangular: unrecognised component name '%s'",
		   compname.c_str());
	return NULL;
    }
	    
    RectangularImperfectPlate *mp = dynamic_cast<RectangularImperfectPlate*>(comp);

    OutputRectangular *op = new OutputRectangular(mp, pan, x, y);    

    logMessage(3, "Creating output %s, %f, %f, %f", compname.c_str(), x, y, pan);
    return op;
}
