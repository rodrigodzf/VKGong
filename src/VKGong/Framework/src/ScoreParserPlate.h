/*
 * VK-Gong code 
 *
 * Author of this version: Àngels Aragonès
 *
 * Based on VK-Gong for Matlab and NeSS Framework Code (Copyright (c) The University of Edinburgh, 2014-2015. All rights reserved. Author: James Perry)
 *
 * Score file parser for the plate format. It has two subclasses for circular and rectangular plates. 
 */
#ifndef _SCORE_PARSER_PLATE_H_
#define _SCORE_PARSER_PLATE_H_

#include "ScoreParser.h"


#include <iostream>
using namespace std;

class ScoreParserPlate : public ScoreParser {
 public:
    ScoreParserPlate(string filename);
    virtual ~ScoreParserPlate();
    double getDuration();

 protected:
    virtual bool parse(PhysicalSystem *physicalSystem);

    Input *parseStrike(istream &in);
    Input *parseHarmonicInput(istream &in);
    Input *parseNoiseInput(istream &in);
 

    PhysicalSystem *physicalSystem;
    bool durationOnly;
    double duration;

    virtual int handleItem(string type, istream &in);
};

#endif
