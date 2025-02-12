/*
 * VK-Gong code 
 *
 * Based on VK-Gong for Matlab and NeSS Framework Code (Copyright (c) The University of Edinburgh, 2014-2015. All rights reserved. Author: James Perry)
 *
 */

#include "WavWriter.h"
#include "GlobalSettings.h"
#include "Logger.h"

#include <iostream>
#include <fstream>
using namespace std;

WavWriter::WavWriter(string filename)
{
    logMessage(1, "Initialising WavWriter for %s", filename.c_str());
    this->filename = filename;
}

WavWriter::~WavWriter()
{
}

bool WavWriter::writeMonoWavFile(Output *output, double max)
{
    logMessage(1, "Writing mono wav file to %s", filename.c_str());
    GlobalSettings *gs = GlobalSettings::getInstance();

    int sr = (int)gs->getSampleRate();
    int len = gs->getNumTimesteps();

    short *buffer = new short[len];

    double *data = output->getData();

    double maxout = (gs->getMaxOut() * 16384.0);
    
    ofstream fff("./OutputResults/dataWAV.bin", ios::out | ios::binary);
	if (fff.good()) {
		fff.write((const char *)data, len*sizeof(double));
		fff.close();
	}

    // normalise level
    for (int i = 0; i < len; i++) {
	buffer[i] = (short)((data[i] / max) * maxout);
    }

    bool result = writeMonoWav(buffer, len, sr);

    delete[] buffer;
    return result;
}

bool WavWriter::writeStereoMix(vector<Output*> *outputs, double max)
{
    logMessage(1, "Writing stereo mix wav file to %s", filename.c_str());
    GlobalSettings *gs = GlobalSettings::getInstance();

    int sr = (int)gs->getSampleRate();
    int len = gs->getNumTimesteps();

    short *buffer = new short[len * 2];

    double maxout = (gs->getMaxOut() * 16384.0);

    // mix and normalise
    for (int i = 0; i < len; i++) {
	double left = 0.0;
	double right = 0.0;
	for (int j = 0; j < outputs->size(); j++) {
	    // mix according to pan for that output
	    left += (1.0 - ((outputs->at(j)->getPan() + 1.0) / 2.0)) *
		outputs->at(j)->getData()[i];
	    right += ((outputs->at(j)->getPan() + 1.0) / 2.0) *
		outputs->at(j)->getData()[i];
	}
	// normalise level
	buffer[i*2] = (short)((left / max) * maxout);
	buffer[i*2+1] = (short)((right / max) * maxout);
    }

    bool result = writeStereoWav(buffer, len, sr);

    delete[] buffer;
    return result;
}

void WavWriter::initWavHeader(WavHeader *hdr, bool stereo, int sampleRate,
			      int len)
{
    logMessage(1, "initialising WavHeader: %d, %d, %d", (int)stereo, sampleRate, len);
    // FIXME: this won't work on a big endian system
    hdr->chunkID[0] = 'R';
    hdr->chunkID[1] = 'I';
    hdr->chunkID[2] = 'F';
    hdr->chunkID[3] = 'F';

    int nChannels = 1;
    if (stereo) nChannels = 2;

    hdr->chunkSize = 36 + (len * nChannels * 2);

    hdr->format[0] = 'W';
    hdr->format[1] = 'A';
    hdr->format[2] = 'V';
    hdr->format[3] = 'E';

    hdr->subChunkID[0] = 'f';
    hdr->subChunkID[1] = 'm';
    hdr->subChunkID[2] = 't';
    hdr->subChunkID[3] = ' ';

    hdr->subChunkSize = 16;
    hdr->audioFormat = 1;
    hdr->numChannels = nChannels;
    hdr->sampleRate = sampleRate;
    hdr->byteRate = sampleRate * nChannels * 2;
    hdr->blockAlign = nChannels * 2;
    hdr->bitsPerSample = 16;

    hdr->subChunkID2[0] = 'd';
    hdr->subChunkID2[1] = 'a';
    hdr->subChunkID2[2] = 't';
    hdr->subChunkID2[3] = 'a';

    hdr->subChunkSize2 = len * nChannels * 2;
}

bool WavWriter::writeMonoWav(short *data, int len, int sr)
{
    WavHeader hdr;

    initWavHeader(&hdr, false, sr, len);

    ofstream of(filename.c_str(), ios::out | ios::binary);
    if (!of.good()) {
	logMessage(5, "Error opening file %s for writing", filename.c_str());
	return false;
    }

    of.write((const char *)&hdr, sizeof(hdr));
    of.write((const char *)data, len * 2);

    of.close();
    return true;
}

bool WavWriter::writeStereoWav(short *data, int len, int sr)
{
    WavHeader hdr;

    initWavHeader(&hdr, true, sr, len);

    ofstream of(filename.c_str(), ios::out | ios::binary);
    if (!of.good()) {
	logMessage(5, "Error opening file %s for writing", filename.c_str());
	return false;
    }

    of.write((const char *)&hdr, sizeof(hdr));
    of.write((const char *)data, len * 2 * 2);

    of.close();
    return true;
}
