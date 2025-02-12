/*
 * VK-Gong code 
 *
 * Author of this version: Àngels Aragonès
 *
 * Based on VK-Gong for Matlab and NeSS Framework Code (Copyright (c) The University of Edinburgh, 2014-2015. All rights reserved. Author: James Perry)
 *
 * This class represents any component within a simulation, with subclasses
 * for circular and rectangular imperfect plates. It assumes that every
 * component has 3 state arrays containing its state at the most recent
 * 3 timesteps. Any component that needs more than this will have to manage
 * them itself in a subclass.
 */
#ifndef _COMPONENT_H_
#define _COMPONENT_H_


extern "C" {
#include "csrmatrix.h"
#include "banded.h"
};


//#include "Task.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
using namespace std;

class Input;


class Component {
 public:
    Component(string name);
    virtual ~Component();

    // getters for state arrays and overall state size
    virtual double *getU() { return u; }
    virtual double *getU1() { return u1; }
    virtual double *getU2() { return u2; }
    int getStateSize() { return Nphi; }

    string getName() { return name; }


    // component hierarchy
    Component *getParent() { return parent; }
    vector<Component*> *getChildren() { return &children; }

    void addInput(Input *input) { inputs.push_back(input); }
    vector<Input*> *getInputs() { return &inputs; }

    // swap the component's state buffers, used at end of timestep
    virtual void swapBuffers(int n);
    virtual void swapBuffersEta(int n);

    // run one timestep (abstract, specific to component type)
    virtual void runTimestepVerlet(int n);
    virtual void runTimestepECS(int n);
 

 protected:
    void runInputs(int n, double *s, double *s1, double *s2);

    void doSaveState();

    // name
    string name;

    // state arrays
    double *u, *u1, *u2;
    
    // Time integration variables
    double *Ai;
	
	double *t0, *t0t, *t1, *t2, *t3, *t4, *t5, *t7, *t8;
	double *H0, *H1, *H2;
	double *qAi, *q2Ai, *qq1;
	double *eta_temp, *eta, *eta1, *eta2;
	double *G, *Ga, *Gb;
	double *G0, *Gc1, *Gc2, *Gc;
	double *mat_imp;
	double *C, *C1, *C2;
    
    int Npsi;
    int Nphi;

    // overall state size
    //int ss;


    // hierarchical relationships
    Component *parent;
    vector<Component*> children;

    // inputs targetting this component
    vector<Input*> inputs;

    // save out state every timestep!
    bool logState;

    // file for logging state to
    ofstream *stateStream;
};

#endif
