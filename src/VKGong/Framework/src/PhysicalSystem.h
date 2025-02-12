/*
 * VK-Gong code 
 *
 * Author of this version: Àngels Aragonès
 *
 * Based on VK-Gong for Matlab and NeSS Framework Code (Copyright (c) The University of Edinburgh, 2014-2015. All rights reserved. Author: James Perry)
 *
 * This class represents an physical system, which is a complete set of components
 * and other entities used to run a simulation.
 */
#ifndef _PHYSICAL_SYSTEM_H_
#define _PHYSICAL_SYSTEM_H_

#include "Component.h"
#include "Output.h"
#include "Profiler.h"

#include <vector>
#include <string>
using namespace std;

class PhysicalSystem {
 public:
    PhysicalSystem();
    virtual ~PhysicalSystem();

    void addComponent(Component *comp) { components.push_back(comp); }
    void addOutput(Output *output) { outputs.push_back(output); }

    virtual void runTimestepVerlet(int n);
    virtual void runTimestepECS(int n);
    virtual void endTimestepVerlet(int n);
    virtual void endTimestepECS(int n);

    void saveOutputs(string outputname, bool individual, bool raw);

    Component *getComponentByName(string name);

    vector<Component*> *getComponents() { return &components; }
    vector<Output*> *getOutputs() { return &outputs; }


 protected:
    double getMaxOutput();

    // all the components
    vector<Component*> components;
    vector<Output*> outputs;

    Profiler *profiler;
};

#endif
