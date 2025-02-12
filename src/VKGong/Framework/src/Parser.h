/*
 * VK-Gong code 
 *
 * Author of this version: Àngels Aragonès
 *
 * Based on VK-Gong for Matlab and NeSS Framework Code (Copyright (c) The University of Edinburgh, 2014-2015. All rights reserved. Author: James Perry)
 *
 * Parser superclass, providing text file parsing abilities to the instrument and score
 * parsers.
 */

#ifndef _PARSER_H_
#define _PARSER_H_

#include <string>
#include <iostream>
using namespace std;

class Parser {
 public:
    Parser(string filename);
    virtual ~Parser();

 protected:
    string filename;

    int parseTextFile();

    virtual int handleItem(string type, istream &in);
};

#endif
