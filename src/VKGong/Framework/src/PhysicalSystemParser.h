/*
 * VK-Gong code 
 *
 * Author of this version: Àngels Aragonès
 *
 * Based on VK-Gong for Matlab and NeSS Framework Code (Copyright (c) The University of Edinburgh, 2014-2015. All rights reserved. Author: James Perry)
 *
 * Abstract superclass for PhysicalSystem file parsers.
 */
#ifndef _PHYSICAL_SYSTEM_PARSER_H_
#define _PHYSICAL_SYSTEM_PARSER_H_

#include "Parser.h"
#include "PhysicalSystem.h"

#include <string>
using namespace std;

class PhysicalSystemParser : public Parser {
 public:
    PhysicalSystemParser(string filename);
    virtual ~PhysicalSystemParser();

    static PhysicalSystem *parsePhysicalSystem(string filename);

 protected:
    PhysicalSystem *physicalSystem;
    virtual PhysicalSystem *parse() = 0;
};

#endif
