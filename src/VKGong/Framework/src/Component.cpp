/*
 * VK-Gong code 
 *
 * Author of this version: Àngels Aragonès
 *
 * Based on VK-Gong for Matlab and NeSS Framework Code (Copyright (c) The University of Edinburgh, 2014-2015. All rights reserved. Author: James Perry)
 *
 */
#include "Component.h"
#include "Input.h"
#include "SettingsManager.h"
#include "GlobalSettings.h"
#include "Logger.h"
#include "MathUtil.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cstring>
using namespace std;

Component::Component(string name)
{
    this->name = name;
    parent = NULL;
    u = NULL;
    u1 = NULL;
    u2 = NULL;
    
    Ai = NULL;
    t0 = NULL;
    t0t = NULL;
    t1= NULL;
	t2 = NULL;
	t3 = NULL;
	t4 = NULL;
	t5 = NULL;
	t7 = NULL;
	t8 = NULL;
	
	H0 = NULL;
	H1 = NULL;
	H2 = NULL;
	
	qAi = NULL;
	q2Ai = NULL;
	qq1 = NULL;
	
	eta_temp = NULL;
	eta = NULL;
	eta1 = NULL;
	eta2 = NULL;
	
	G = NULL;
	Ga = NULL;
	Gb = NULL;
	
	G0 = NULL;
	Gc1 = NULL;
	Gc2 = NULL;
	Gc = NULL;
	
	mat_imp = NULL;
	
	C = NULL;
	C1 = NULL;
	C2 = NULL;
	
	

    stateStream = NULL;
    logState = SettingsManager::getInstance()->getBoolSetting(name, "log_state");
}

Component::~Component()
{

    int i;
    if (u) delete[] u;
    if (u1) delete[] u1;
    if (u2) delete[] u2;
    if (Ai) delete[] Ai;
    if (t0) delete[] t0;
    if (t0t) delete[] t0t;
    if (t1) delete[] t1;
    if (t2) delete[] t2;
    if (t3) delete[] t3;

    if (t4) delete[] t4;
    if (t5) delete[] t5;
    if (t7) delete[] t7;
    if (t8) delete[] t8;

    if (H0) delete[] H0;
    if (H1) delete[] H1;
    if (H2) delete[] H2;
	
    if (qAi) delete[] qAi;
    if (q2Ai) delete[] q2Ai;
    if (qq1) delete[] qq1;
  
    if (eta_temp) delete[] eta_temp;
    if (eta2) delete[] eta2;
    if (eta1) delete[] eta1;
    if (eta) delete[] eta;
			
    if (G) delete[] G;
    if (Ga) delete[] Ga;
    if (Gb) delete[] Gb;
	
    if (G0) delete[] G0;
    if (Gc1) delete[] Gc1;
    if (Gc2) delete[] Gc2;
    if (Gc) delete[] Gc;
	
    if (mat_imp) delete[] mat_imp;
	
    if (C) delete[] C;
    if (C1) delete[] C1;
    if (C2) delete[] C2;
    
    for (i = 0; i < inputs.size(); i++) {
	delete inputs.at(i);
    }

    if (stateStream != NULL) {
	stateStream->close();
	delete stateStream;
    }
    


}

void Component::doSaveState()
{
    if (logState) {
	if (stateStream == NULL) {
	    char filenamebuf[1000];
	    sprintf(filenamebuf, "%s-state.bin", name.c_str());
	    stateStream = new ofstream(filenamebuf, ios::out | ios::binary);
	}
	stateStream->write((const char *)u, Nphi * sizeof(double));
    }
}

void Component::swapBuffers(int n)
{
    double *tmp;

    doSaveState();

    tmp = u2;
    u2 = u1;
    u1 = u;
    u = tmp;
}

void Component::swapBuffersEta(int n)
{
    double *tmp;

    tmp = eta2;
    eta2 = eta1;
    eta1 = eta;
    eta = tmp;
}

void Component::runInputs(int n, double *s, double *s1, double *s2)
{
    int i;
    
    for (i = 0; i < inputs.size(); i++) {

	    inputs.at(i)->runTimestep(n, s, s1, s2);

    }
}

void Component::runTimestepVerlet(int n){
    int i;

	// run standard version
	/* t0 = (H1*(q1+Ai)); */
    for (i = 0; i < Nphi; i++) {
    	qAi[i] = u1[i] + Ai[i]; 	 

    }
    
    
	denseMatrixVectorMultiply(H1, qAi, t0, (Npsi*Nphi), Nphi);
	

    /* t1 = (H1*q1); */
	denseMatrixVectorMultiply(H1, u1, t1, (Npsi*Nphi), Nphi);
	
    /* t2 = t1*(q1 + 2Ai); */
    for (i = 0; i < Nphi; i++) {
    	q2Ai[i] = u1[i] + 2*Ai[i]; 	 	
    }
    denseMatrixVectorMultiply(t1, q2Ai, t2, Npsi, Nphi);
    
    /* G = t0.'*t2; */
    denseMatrixVectorMultiplyTransposed(t0, t2, G, Npsi, Nphi);
    
    for (i = 0; i < Nphi; i++) {
    	G[i] = G[i] / C[i];
    }

    /* q = - C1.*q1 - C2.*q2 - G + inputs; */
    for (i = 0; i < Nphi; i++) {
	u[i] = -(C1[i]*u1[i]) - (C2[i]*u2[i]) - G[i];

    }
    
    runInputs(n, u, u1, u2);
    

    

}

void Component::runTimestepECS(int n){
    int i,j;
   

	// t0 = (H1*(q1+Ai)); 
    for (i = 0; i < Nphi; i++) {
    	qAi[i] = u1[i] + Ai[i]; 
    }
	denseMatrixVectorMultiply(H1, qAi, t0, (Npsi*Nphi), Nphi);

	
	transposeMatrixRectangular(Npsi, Nphi, t0, t0t);
	
	// G0 = t0.'*t0 
    denseMatrixMatrixMultiplyRectangular(Nphi, Npsi, t0, t0t, G0); // G0 is symmetric
 	
	// mat_imp = C + G0 
	for (i = 0; i < Nphi; i++) {
		for (j = 0; j < Nphi; j++) {
	    	if (i==j){
	    		mat_imp[i + j*Nphi] = C[i] + G0[i + j*Nphi]; 
	    	}
	    	else{
	    		mat_imp[i + j*Nphi] = G0[i + j*Nphi];
	    	}

		}
	}
	

	
	 //t1 = (H1*q1); 
	denseMatrixVectorMultiply(H1, u1, t1, (Npsi*Nphi), Nphi);
	


	// t4 = t1*Ai; 
	denseMatrixVectorMultiply(t1, Ai, t4, Npsi, Nphi);

	// Ga = t0.'*t4; 
	denseMatrixVectorMultiplyTransposed(t0, t4, Ga, Npsi, Nphi);

	// eta_temp = eta2 - eta1;
	for (i = 0; i < Npsi; i++) {
		eta_temp[i] = eta2[i] - eta1[i];

	}

	 //t5 = (H0*(q1 + Ai)); 
	denseMatrixVectorMultiply(H0, qAi, t5, (Npsi*Nphi), Nphi);

	 
	//Gb = (t5.'*eta_temp); 
	denseMatrixVectorMultiplyTransposed(t5, eta_temp, Gb, Npsi, Nphi); // Repassar que estigui be
	
	
	//G = -C1.*q1 - C2.*q2 - (Ga' - Gb) 
	for (i = 0; i < Nphi; i++) {
		G[i] = -(C1[i]*u1[i]) - (C2[i]*u2[i]) - (Ga[i] - Gb[i]);


		
	}
	

	    
    runInputs(n, G, u1, u2);
    
     
    //q = mat_imp\(- C1.*q1 - C2.*q2 - Gtotal + f_time(:,i));
 
    //if (gaussSolver(mat_imp, G, &minPivot, &row, Nphi)){ // Alerta que treballa per files. tenir-ho en compte en els passos previs. Tamb'e tenir en compte Cnorm. 
    if (gaussSolver(mat_imp, G, Nphi)){    
    	for (i = 0; i < Nphi; i++) {
  			u[i] = G[i]; 
		}
    	
    	// t7 = H2*q;
    	denseMatrixVectorMultiply(H2, u, t7, (Npsi*Nphi), Nphi);

    	
    	 //Gc_1 = t7*q1;
    	denseMatrixVectorMultiply(t7, u1, Gc1, Npsi, Nphi);

    	
    	//t8 = H2*Ai;
    	denseMatrixVectorMultiply(H2, Ai, t8, (Npsi*Nphi), Nphi);

    	
    	// Gc_2 = t8*(q + q1);
    	for (i = 0; i < Nphi; i++) {
  			qq1[i] = u[i] + u1[i]; 
		}
    	denseMatrixVectorMultiply(t8, qq1, Gc2, Npsi, Nphi);

    	    
    	// Gc = (Gc_1 + Gc_2);
    	// eta = -eta1 - Gc;
    	
    	for ( i = 0; i < Npsi; i++) {
  			eta[i] = -eta1[i] - (Gc1[i] + Gc2[i]); 
		}


    	    
    }
    else{
    	logMessage(5, "Error: Matrix could not be inverted.");

    }
    
    
    
    
    
}



