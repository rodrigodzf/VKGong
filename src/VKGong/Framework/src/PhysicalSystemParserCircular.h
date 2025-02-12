/*
 * VK-Gong code 
 *
 * Author of this version: Àngels Aragonès
 *
 * Based on VK-Gong for Matlab and NeSS Framework Code (Copyright (c) The University of Edinburgh, 2014-2015. All rights reserved. Author: James Perry)
 *
 * Parser for the circular plate physical system.
 */
#ifndef _PHYSICAL_SYSTEM_PARSER_CIRCULAR_H_
#define _PHYSICAL_SYSTEM_PARSER_CIRCULAR_H_

#include "PhysicalSystemParser.h"
#include "CircularImperfectPlate.h"
#include "OutputCircular.h"

#include <iostream>
using namespace std;

class PhysicalSystemParserCircular : public PhysicalSystemParser {
 public:
    PhysicalSystemParserCircular(string filename);
    virtual ~PhysicalSystemParserCircular();

 protected:
    virtual PhysicalSystem *parse();
    virtual int handleItem(string type, istream &in);

    CircularImperfectPlate *readCircularImperfectPlate(istream &in);
    OutputCircular *readOutput(istream &in);
};

#endif
