/*
 * VK-Gong code 
 *
 * Based on VK-Gong for Matlab and NeSS Framework Code (Copyright (c) The University of Edinburgh, 2014-2015. All rights reserved. Author: James Perry)
 *
 * Reads WAV file inputs
 */

#ifndef _WAVREADER_H_
#define _WAVREADER_H_

#include <string>
using namespace std;

class WavReader {
 public:
    WavReader(string filename);
    ~WavReader();

    double *getValues();
    int getSize();
    int getSampleRate();

    void resampleTo(int sr);
    void normalise(double gain);

 private:
    double *data;
    int size;
    int sampleRate;
};

#endif
