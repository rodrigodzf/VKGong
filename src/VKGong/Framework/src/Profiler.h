/*
 * VK-Gong code 
 *
 * Author of this version: Àngels Aragonès
 *
 * Based on VK-Gong for Matlab and NeSS Framework Code (Copyright (c) The University of Edinburgh, 2014-2015. All rights reserved. Author: James Perry)
 *
 * Profiler object
 */

#ifndef _PROFILER_H_
#define _PROFILER_H_

#include <sys/time.h>

#include <string>
using namespace std;

class Profiler {
 public:
    Profiler(int size, char **names = NULL);
    Profiler();
    virtual ~Profiler();

    void start(int idx);
    void end(int idx);
    double get(int idx);
    string print();

    double getTime();

 protected:
    void init(int size, char **names);

    char **names;

    double *times;
    double startTime;
    int size;
    int maxUsed;

    struct timeval tp;
};

#endif
