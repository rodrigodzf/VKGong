/*
 * VK-Gong code 
 *
 * Author of this version: Àngels Aragonès
 *
 * Based on VK-Gong for Matlab and NeSS Framework Code (Copyright (c) The University of Edinburgh, 2014-2015. All rights reserved. Author: James Perry)
 *
 * Parser for the imperfect rectangular plate physical system.
 */
#ifndef _PHYSICAL_SYSTEM_PARSER_RECTANGULAR_H_
#define _PHYSICAL_SYSTEM_PARSER_RECTANGULAR_H_

#include "PhysicalSystemParser.h"
#include "RectangularImperfectPlate.h"
#include "OutputRectangular.h"

#include <iostream>
using namespace std;

class PhysicalSystemParserRectangular : public PhysicalSystemParser {
 public:
    PhysicalSystemParserRectangular(string filename);
    virtual ~PhysicalSystemParserRectangular();

 protected:
    virtual PhysicalSystem *parse();
    virtual int handleItem(string type, istream &in);

    RectangularImperfectPlate *readRectangularImperfectPlate(istream &in);
    OutputRectangular *readOutput(istream &in);
};

#endif
