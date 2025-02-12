/*
 * VK-Gong code 
 *
 * Author of this version: Àngels Aragonès
 *
 * Based on VK-Gong for Matlab and NeSS Framework Code (Copyright (c) The University of Edinburgh, 2014-2015. All rights reserved. Author: James Perry)
 *
 */
#include "PhysicalSystemParser.h"



#include "PhysicalSystemParserCircular.h"
#include "PhysicalSystemParserRectangular.h"


#include "Logger.h"

PhysicalSystemParser::PhysicalSystemParser(string filename) : Parser(filename)
{
}

PhysicalSystemParser::~PhysicalSystemParser()
{
}

PhysicalSystem *PhysicalSystemParser::parsePhysicalSystem(string filename)
{
    PhysicalSystem *result = NULL;
    PhysicalSystemParser *parser;
    



    parser = new PhysicalSystemParserCircular(filename);
    result = parser->parse();
    delete parser;
    if (result) return result;

    parser = new PhysicalSystemParserRectangular(filename);
    result = parser->parse();
    delete parser;
    if (result) return result; 
    
    return result;
    

}
