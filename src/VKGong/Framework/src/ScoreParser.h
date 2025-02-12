/*
 * VK-Gong code 
 *
 * Author of this version: Àngels Aragonès
 *
 * Based on VK-Gong for Matlab and NeSS Framework Code (Copyright (c) The University of Edinburgh, 2014-2015. All rights reserved. Author: James Perry)
 *
 * Abstract superclass for score file parsers.
 */
#ifndef _SCORE_PARSER_H_
#define _SCORE_PARSER_H_

#include "Parser.h"
#include "PhysicalSystem.h"

#include <string>
using namespace std;

class ScoreParser : public Parser {
 public:
    virtual ~ScoreParser();

    static bool parseScore(string filename, PhysicalSystem *physicalSystem);

    static double getDuration(string filename);

 protected:
    ScoreParser(string filename);

    virtual bool parse(PhysicalSystem *physicalSystem) = 0;
};

#endif
