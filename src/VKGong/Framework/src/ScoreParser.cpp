/*
 * VK-Gong code 
 *
 * Author of this version: Àngels Aragonès
 *
 * Based on VK-Gong for Matlab and NeSS Framework Code (Copyright (c) The University of Edinburgh, 2014-2015. All rights reserved. Author: James Perry)
 *
 */

#include "ScoreParser.h"



#include "ScoreParserPlate.h"



ScoreParser::ScoreParser(string filename) : Parser(filename)
{
}

ScoreParser::~ScoreParser()
{
}

bool ScoreParser::parseScore(string filename, PhysicalSystem *physicalSystem)
{
    bool result = false;
    ScoreParser *parser;



    parser = new ScoreParserPlate(filename);
    result = parser->parse(physicalSystem);
    delete parser;
    if (result) return result;



    return result;
}

double ScoreParser::getDuration(string filename)
{
    double duration;

   
    
    ScoreParserPlate *parserPlate = new ScoreParserPlate(filename);
    duration = parserPlate->getDuration();
    delete parserPlate;
    if (duration > 0.0) return duration;






    return duration;
}
