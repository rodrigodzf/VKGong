/*
 * VK-Gong code 
 *
 * Author of this version: Àngels Aragonès
 *
 * Based on VK-Gong for Matlab and NeSS Framework Code (Copyright (c) The University of Edinburgh, 2014-2015. All rights reserved. Author: James Perry)
 *
 */

#include "CircularImperfectPlate.h"
#include "Logger.h"
#include "GlobalSettings.h"
#include "MathUtil.h"
#include "SettingsManager.h"
#include "Input.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>

#include <algorithm>

#include <gsl/gsl_sf_bessel.h>

using namespace std;

void getCacheDirectory(char *buf, int bufsize);

CircularImperfectPlate::CircularImperfectPlate(string name, int iNphi, int iNpsi, double Rd, double h, double H, char ImperfectionType, double tau2, int modeType,
	       double nu, double Young, double rho, char BC, double KR, double KT, int Nr, int Nth, double dFac, double dExp, double dCons, int fs)
    : Component(name)
{
    int i;
    GlobalSettings *gs = GlobalSettings::getInstance();

    logMessage(1, "Creating CircularImperfectPlate: modes=%d,%d, Radius=%f, h=%f, Height=%f, nu=%f, E=%f, rho=%f", iNphi, iNpsi, Rd, h, H, nu, Young, rho);

    this->Rd = Rd;
    this->H = H;
    
    this->dFac = dFac;
    this->dExp = dExp;
    this->dCons = dCons;

    Nphi = iNphi;
    Npsi = iNpsi;
    this->nu = nu;
    this->Young = Young;
    this->rho = rho;
    this->h = h;
    this->tau2 = tau2;
    this->ImperfectionType = ImperfectionType;
    this->modeType = modeType;
    this->BC = BC;
    this->KR = KR;
    this->KT = KT;
    if (BC == 'c'){
    	this->KR = -1;
    	this->KT = -1;
    } 
    else if (BC == 'f'){
    	this->KR = 0;
    	this->KT = 0;
    }
    
    
    this->Nr = Nr;
    this->Nth = Nth;
    
    
    this->e = 12*(1-nu*nu);
    D = Young*h*h*h/e;
    
    tnd = Rd*Rd*sqrt(rho*h/D);
    this->tnd = tnd;
 

    // allocate state arrays
    u = new double[Nphi];
    u1 = new double[Nphi];
    u2 = new double[Nphi];
    

    memset(u, 0, Nphi * sizeof(double));
    memset(u1, 0, Nphi * sizeof(double));
    memset(u2, 0, Nphi * sizeof(double));
    
    
    // allocate temporary buffers corresponding to every scheme.
    if (gs->getStormerVerletScheme()){
    	
		t0 = new double[Npsi*Nphi];
		t1 = new double[Npsi*Nphi];
		t2 = new double[Npsi];
		t3 = new double[Npsi];
		qAi = new double[Nphi];
		q2Ai = new double[Nphi];
		G = new double[Nphi];
		
		H1 = new double[Nphi*Nphi*Npsi];
		
		memset(G, 0, Nphi * sizeof(double));

    }
    else{
    	t0 = new double[Npsi*Nphi];
    	t0t = new double[Npsi*Nphi];
		t1 = new double[Npsi*Nphi];
		t4 = new double[Npsi];
		t5 = new double[Npsi*Nphi];
		t7 = new double[Npsi*Nphi];
		t8 = new double[Npsi*Nphi];
		qAi = new double[Nphi];
		qq1 = new double[Nphi];
		
		H0 = new double[Nphi*Nphi*Npsi];
		H1 = new double[Nphi*Nphi*Npsi];
		H2 = new double[Nphi*Nphi*Npsi];
		
		eta = new double[Npsi];
		eta1 = new double[Npsi];
		eta2 = new double[Npsi];
		eta_temp = new double[Npsi];
		
		memset(eta, 0, Npsi * sizeof(double));
		memset(eta1, 0, Npsi * sizeof(double));
		memset(eta2, 0, Npsi * sizeof(double));
		memset(eta_temp, 0, Npsi * sizeof(double));
		
		Ga = new double[Nphi];
		Gb = new double[Nphi];
		G = new double[Nphi];
		
		memset(G, 0, Nphi * sizeof(double));
		memset(Ga, 0, Nphi * sizeof(double));
		memset(Gb, 0, Nphi * sizeof(double));
		
		Gc1 = new double[Npsi];
		Gc2 = new double[Npsi];
		Gc = new double[Npsi];
		
		memset(Gc, 0, Npsi * sizeof(double));
		memset(Gc1, 0, Npsi * sizeof(double));
		memset(Gc2, 0, Npsi * sizeof(double));
		
		G0 = new double[Nphi*Nphi];
		mat_imp = new double[Nphi*Nphi];
		
		memset(G0, 0, Nphi * Nphi * sizeof(double));
		memset(mat_imp, 0, Nphi * Nphi * sizeof(double));
		
    }
  

    
   
    double dx=1e-3;
    double xmax=100;
    double dr_H=1e-4;
    
    xi = new double[Nphi];
	kt = new int[Nphi];
	ct = new int[Nphi];
	ov = new double[Nphi];
	zl = new double[Npsi];
	kl = new int[Npsi];
	cl = new int[Npsi];

	
    TransverseModeCalculator(dx, xmax, BC, nu, KR, KT, Nphi, xi, kt, ct, ov); 
  
    InPlaneModeCalculator(dx, xmax, BC, nu, Npsi, zl, kl, cl);
    
    LoadHTensor(dr_H);
    
	
    
    Ai = new double[Nphi];
    memset(Ai, 0, Nphi * sizeof(double));
    
    if (H>0){

    // see if cached Ai for this imperfection exists
    char cachefile[1000];
       sprintf(cachefile, "./CircularPlateData/Imperfection_Ai_%d_%f_%f_%c_%c_%f.bin", Nphi, H,Rd,ImperfectionType, BC, nu); 
       ifstream fi2(cachefile, ios::in | ios::binary);
       if (fi2.good()) {
    	   logMessage(1, "Loading Ai from cache");
   	// found cache file, load it
    	   fi2.read((char *)Ai, Nphi*sizeof(double));
    	   fi2.close();
       }
       else {
    	   double *Imperfection = new double[Nr*Nth*sizeof(double)];  
    	   
    	   AxisymmetricCap(H,Rd, ImperfectionType, Imperfection, Nr, Nth);    	   
    	   double max_error = 0.1;
  	   
    	   ProjectionCoefficients(Imperfection, Ai, Nphi, modeType, max_error, Nr, Nth, nu, KR, kt, nt, ct, xi);
    	      	   
    	   logMessage(1, "Performed ProjectionCoefficients");

   	// write to cache for next time
    	   ofstream fo2(cachefile, ios::out | ios::binary);
    	   if (fo2.good()) {
    		   fo2.write((const char *)Ai, Nphi*sizeof(double));
    		   fo2.close();
    	   }
    	   delete[] Imperfection;
       }
    }

	
    C = new double[Nphi];
    C1 = new double[Nphi];
    C2 = new double[Nphi];
    Cnorm = new double[Nphi];
       
    
    plate_def(&fs);

    logMessage(1, "Performed plate_def, sample rate=%d", fs);

    
    gs->setSampleRate((double)fs);

    logMessage(1, "Final mode counts: A=%d, Nphi=%d", this->Npsi, this->Nphi);
  
}

CircularImperfectPlate::~CircularImperfectPlate()
{

    delete[] Cnorm;
    delete[] ov;
    delete[] xi;
    delete[] kt;
    delete[] ct;	
    delete[] zl;
    delete[] kl;
    delete[] cl;

}


void CircularImperfectPlate::plate_def(int *fs)
{
    double fs_lim;
    double fsd;
    int m;
    double *c = new double[Nphi];
    double c_coef;
    
    GlobalSettings *gs = GlobalSettings::getInstance();
   
    fs_lim = ov[Nphi-1] / tnd / 2.0; // ov is already dimensionless fs = fsd * tnd; ov_d = ov / tnd;
	if (*fs < fs_lim) {    
		*fs = (int)(fs_lim * 3); // fs is stored in dimensional form
		printf("Warning: The sampling rate introduced by user is too low. It will be modified to fulfill the stability limit. \n");
	}

    printf("Simulation fs = %d \n", *fs);
    
    fsd = (double)*fs; // fsd stands for fs double
    fsd = round(fsd*tnd);

    /* damping ratios */
    
    c_coef = Rd*Rd*rho*h/sqrt(rho*h*D);
    
    for (m = 0; m < Nphi; m++) {

    	c[m] = dFac * pow(ov[m]/tnd, dExp) + dCons;
    	c[m] = c_coef*c[m];
    	
    }

    for (m = 0; m < Nphi; m++) {
	C[m] = (fsd*fsd + c[m]*fsd); 
	C1[m] = (-2.0*fsd*fsd + ov[m]*ov[m]);
	C2[m] = (fsd*fsd - c[m]*fsd);

	
	if (gs->getStormerVerletScheme()){
			Cnorm[m] = C[m]; /* Cnorm is just a copy of C used if StormerVerlet */
			C1[m] = C1[m] / C[m];
			C2[m] = C2[m] / C[m];
	}
	else{
		Cnorm[m] = 1.0;
	}
		
    }
    delete[] c;
}



void CircularImperfectPlate::TransverseModeCalculator(double dx, double xmax, char BC, double nu, double KR, double KT, int Nphi, double* xi, int* kt, int* ct, double* ov )
{
	// Transverse eigenfrequencies
    double *mode_t;
    char filename[1000];
    int length, width = 6;


    switch (BC){
    	case 'f':
    		sprintf(filename, "./CircularPlateData/mode_t_free-nu_%f.bin",nu);
    		break;
    	case 'e':
    		sprintf(filename, "./CircularPlateData/mode_t_elastic-nu_%f-KR_%f-KT_%f.bin", nu, KR, KT);
    		break;
    	case 'c':
    		sprintf(filename, "./CircularPlateData/mode_t_clamped.bin");
    		break;
    }
    
	ifstream fi(filename, ios::in | ios::binary);

	

    if (fi.good()) { // found file, load it

    	logMessage(1, "Loading transverse eigenfrequencies file");

        fi.seekg (0, fi.end);
        length = fi.tellg();
        fi.seekg (0, fi.beg);

        length = length/width/8;


    	mode_t = new double[length*width];

    	fi.read((char *)mode_t, length*width*sizeof(double));
    	fi.close();

        for (int n=0; n<Nphi; n++){
        	        	
        	xi[n] = mode_t[n*width+1];
        	kt[n] = (int)mode_t[n*width+2];
        	//nt[n] = (int)mode_t[n*width+3];
        	ct[n] = (int)mode_t[n*width+4];
        	ov[n] = mode_t[n*width+5];


        }

        delete[] mode_t;

    }
    else{
    	bool FoundZero = true;
    	int k = 0, n;
    	int j = 0, NumFoundModes=0;
    	double Nx;
    	int  xx_min=0;
    	double x = 0, J3, J2, J1, J0, I3, I2, I1, I0, Jtild, Itild, Jtild2, Itild2;
    	double f, f_prev = 0;
    	int auxArrayLength = 2000;

    	std::vector<CircularImperfectPlate::ModeXKN> zeros(auxArrayLength);
    	
    	Nx = round(xmax/dx);
    	
    	while (FoundZero){
    		FoundZero = false;
    		f_prev = 0;
    		n = 1;
			for (int xx = xx_min; xx<=Nx; xx++){
				x = xx*dx;
				
				if ((BC == 'e')||(BC == 'f' )){
					
    				J3 = gsl_sf_bessel_Jn(k-3,x);   				
					J2 = gsl_sf_bessel_Jn(k-2,x);					
					I3 = gsl_sf_bessel_In(k-3,x);
					I2 = gsl_sf_bessel_In(k-2,x);
					
					if ((fabs(2*(k-2)*J2/x - J3) < (1e-100)) || (fabs(-2*(k-2)*I2/x + I3) < (1e-100))) { // To avoid underflow.
						
						xx_min = xx + 1;
						J3 = 0;
						J2 = 0;
						J1 = 0;
						J0 = 0;
						I3 = 0;
						I2 = 0;
						I1 = 0;
						I0 = 0;
					}
					else {
						J1 = gsl_sf_bessel_Jn(k-1,x);
						J0 = gsl_sf_bessel_Jn(k,x);
						I1 = gsl_sf_bessel_In(k-1,x);
						I0 = gsl_sf_bessel_In(k,x);
					}
					

					Jtild = x*x*J2 + ((nu-2*k+1) + KR)*x*J1 + (k*(k+1)*(1-nu) - KR*k)*J0;
					Itild = x*x*I2 + ((nu-2*k+1) + KR)*x*I1 + (k*(k+1)*(1-nu) - KR*k)*I0;
					
					Jtild2 = x*x*x*J3 + (4-3*k)*x*x*J2 + k*(k*(1+nu)-2)*x*J1 + ((1-nu)*k*k*(1+k) - KT)*J0;
					Itild2 = x*x*x*I3 + (4-3*k)*x*x*I2 + k*(k*(1+nu)-2)*x*I1 + ((1-nu)*k*k*(1+k) - KT)*I0;
					
									
					f = Jtild*Itild2 - Itild*Jtild2;
				}else{
					
					J1 = gsl_sf_bessel_Jn(k-1,x);
					J0 = gsl_sf_bessel_Jn(k,x);
					
					I1 = gsl_sf_bessel_In(k-1,x);
					I0 = gsl_sf_bessel_In(k,x);
					
					if((fabs(2*k*J0/x - J1) < (1e-100)) || (fabs(-2*k*I0/x + I1) < (1e-100))){ // To avoid underflow.
						xx_min = xx + 1;
					}
					
					f = J1*I0 - J0*I1;
					

					
				}
				
							
				if ((f_prev*f)<0){
					FoundZero = true;
					zeros[j].x = ((x-dx)*f - x*f_prev)/(f - f_prev);
					zeros[j].k = k;
					zeros[j].n = n;
 					n++;
					j++;	
					NumFoundModes++;
					if (k>0){NumFoundModes++;}					
				}
				
				f_prev = f;

    		}
			
					
			k++; 	
			
			
    			    		
    	}
    	   	
    	mode_t = SortZeros(j, NumFoundModes, zeros, BC);
    	GlobalSettings *gs = GlobalSettings::getInstance();
    	
    	if (NumFoundModes<Nphi){
    		logMessage(5, "The number of found modes is less than Nphi. Please, increase xmax.");
    		exit(1);
    	}
    	
		for (int n=0; n<Nphi; n++){
		        	        	
			xi[n] = mode_t[n*width+1];
			kt[n] = (int)mode_t[n*width+2];
			//nt[n] = (int)mode_t[n*width+3];
			ct[n] = (int)mode_t[n*width+4];
			ov[n] = mode_t[n*width+5];
		}
		
		logMessage(1, "Performed Transverse mode calculator");
	    
		ofstream fo(filename, ios::out | ios::binary);
		if (fo.good()) {
		    fo.write((const char *)mode_t, NumFoundModes*width*sizeof(double));
		    fo.close();
		}
		
    }

}

void CircularImperfectPlate::InPlaneModeCalculator(double dx, double xmax, char BC, double nu, int Npsi, double* zl, int* kl, int* cl)
{
	// In Plane eigenfrequencies
    double *mode_l;
    char filename[1000];
    int length, width = 6;
    
    
    switch (BC){
    	case 'f':
    		sprintf(filename, "./CircularPlateData/mode_l_free.bin");
    		break;
    	case 'e':
    		sprintf(filename, "./CircularPlateData/mode_l_elastic.bin");
    		break;
    	case 'c':
    		sprintf(filename, "./CircularPlateData/mode_l_clamped-nu_%f.bin", nu);
    		break;
    }
    
    
	ifstream fi(filename, ios::in | ios::binary);


    if (fi.good()) { // found file, load it

    	logMessage(1, "Loading in-plane eigenfrequencies file");

        fi.seekg (0, fi.end);
        length = fi.tellg();
        fi.seekg (0, fi.beg);

        length = length/width/8;

    	mode_l = new double[length*width];

    	fi.read((char *)mode_l, length*width*sizeof(double));
    	fi.close();

        for (int n=0; n<Npsi; n++){
        	zl[n] = mode_l[n*width+1];
        	kl[n] = (int)mode_l[n*width+2];
        	//nl[n] = (int)mode_l[n*width+3];
        	cl[n] = (int)mode_l[n*width+4];

        }

        delete[] mode_l;

    }
    else{
    	bool FoundZero = true;
    	int k = 0, n;
    	int j = 0, NumFoundModes=0;
    	double Nx;
    	int xx_min = 0;
    	double x = 0, J3, J2, J1, J0, I3, I2, I1, I0, Jtild, Itild;
    	double f, f_prev = 0;
    	int auxArrayLength = 2000;

    	std::vector<CircularImperfectPlate::ModeXKN> zeros(auxArrayLength);
    	
    	
    	Nx = round(xmax/dx);
    	
    	while (FoundZero){
    		FoundZero = false;
    		f_prev = 0;
    		n = 1;
			for (int xx = xx_min; xx<=Nx; xx++){
  				x = xx*dx;
				if (BC == 'c'){
      				
					J3 = gsl_sf_bessel_Jn(k-3,x);   				
					J2 = gsl_sf_bessel_Jn(k-2,x);					
					I3 = gsl_sf_bessel_In(k-3,x);
					I2 = gsl_sf_bessel_In(k-2,x);
					
					if ((fabs(2*(k-2)*J2/x - J3) < (1e-100)) || (fabs(-2*(k-2)*I2/x + I3) < (1e-100))) { // To avoid underflow.
						
						xx_min = xx + 1;
						J3 = 0;
						J2 = 0;
						J1 = 0;
						J0 = 0;
						I3 = 0;
						I2 = 0;
						I1 = 0;
						I0 = 0;
					}
					else {
						J1 = gsl_sf_bessel_Jn(k-1,x);
						J0 = gsl_sf_bessel_Jn(k,x);
						I1 = gsl_sf_bessel_In(k-1,x);
						I0 = gsl_sf_bessel_In(k,x);
					}
					
					Jtild = x*x*J2 + (-nu-2*k+1)*x*J1 + (k*(k+1) + nu*k*(1-k))*J0;
					Itild = x*x*I2 + (-nu-2*k+1)*x*I1 + (k*(k+1) + nu*k*(1-k))*I0;
					
					
					 f=Itild*(x*x*x*J3+(4-3*k)*x*x*J2+k*(k*(1-nu)-2)*x*J1+(1+nu)*k*k*(1+k)*J0)-Jtild*(x*x*x*I3+(4-3*k)*x*x*I2+k*(k*(1-nu)-2)*x*I1+(1+nu)*k*k*(1+k)*I0);
					 
				}else{
					
					J1 = gsl_sf_bessel_Jn(k-1,x);
					J0 = gsl_sf_bessel_Jn(k,x);
					
					I1 = gsl_sf_bessel_In(k-1,x);
					I0 = gsl_sf_bessel_In(k,x);
					
					if((fabs(2*k*J0/x - J1) < (1e-100)) || (fabs(-2*k*I0/x + I1) < (1e-100))){ // To avoid underflow.
						xx_min = xx + 1;
					}
					
					f = J1*I0 - J0*I1;									
				}
							
				if ((f_prev*f)<0){
					FoundZero = true;
					zeros[j].x = ((x-dx)*f - x*f_prev)/(f - f_prev);
					zeros[j].k = k;
					zeros[j].n = n;
 					n++;
					j++;	
					
					NumFoundModes++;
					if (k>0){NumFoundModes++;}					
				}
				
				f_prev = f;

    		}
			k++; 			
    			    		
    	}
    	
    	mode_l = SortZeros(j, NumFoundModes, zeros, BC);
    	
    	if (NumFoundModes<Npsi){
    		logMessage(5, "The number of found modes is less than Npsi. Please, increase xmax.");
    		exit(1);
    	}
    	
		for (int n=0; n<Npsi; n++){		        	        	
			zl[n] = mode_l[n*width+1];
			kl[n] = (int)mode_l[n*width+2];
			//nl[n] = (int)mode_l[n*width+3];
			cl[n] = (int)mode_l[n*width+4];

		}
		
		logMessage(1, "Performed InPlane mode calculator");
		// write to cache for next time
	    
		ofstream fo(filename, ios::out | ios::binary);
		if (fo.good()) {
		    fo.write((const char *)mode_l, NumFoundModes*width*sizeof(double));
		    fo.close();
		}
    	

    }

}


#define PI2 (M_PI*M_PI)
#define PI3 (M_PI*M_PI*M_PI)
#define PI4 (M_PI*M_PI*M_PI*M_PI)
#define PI5 (M_PI*M_PI*M_PI*M_PI*M_PI)

void CircularImperfectPlate::LoadHTensor(double dr_H) 
{
	char filenameH0[1000], filenameH1[1000], filenameH2[1000];		
	int NphiF, NpsiF;
	bool HLoaded = false;
		
	GlobalSettings *gs = GlobalSettings::getInstance();
	
	NphiF = gs->getNphiMin();
	NpsiF = gs->getNpsiMin();	
	
	if (Nphi>NphiF){ NphiF = Nphi;}
	if (Npsi>NpsiF){ NpsiF = Npsi;}
	
	if (BC == 'f'){
		sprintf(filenameH0, "./CircularPlateData/H0_free-Nphi_%d-Npsi_%d-nu_%f-dr_%f.bin", NphiF, NpsiF, nu, dr_H);
		sprintf(filenameH1, "./CircularPlateData/H1_free-Nphi_%d-Npsi_%d-nu_%f-dr_%f.bin", NphiF, NpsiF, nu, dr_H);
		sprintf(filenameH2, "./CircularPlateData/H2_free-Nphi_%d-Npsi_%d-nu_%f-dr_%f.bin", NphiF, NpsiF, nu, dr_H);
	} else if (BC == 'c'){
		sprintf(filenameH0, "./CircularPlateData/H0_clamped-Nphi_%d-Npsi_%d-nu_%f-dr_%f.bin", NphiF, NpsiF, nu, dr_H);
		sprintf(filenameH1, "./CircularPlateData/H1_clamped-Nphi_%d-Npsi_%d-nu_%f-dr_%f.bin", NphiF, NpsiF, nu, dr_H);
		sprintf(filenameH2, "./CircularPlateData/H2_clamped-Nphi_%d-Npsi_%d-nu_%f-dr_%f.bin", NphiF, NpsiF, nu, dr_H);
	} else {
		sprintf(filenameH0, "./CircularPlateData/H0_elastic-Nphi_%d-Npsi_%d-nu_%f-dr_%f-KR_%f-KT_%f.bin", NphiF, NpsiF, nu, dr_H, KR, KT);
		sprintf(filenameH1, "./CircularPlateData/H1_elastic-Nphi_%d-Npsi_%d-nu_%f-dr_%f-KR_%f-KT_%f.bin", NphiF, NpsiF, nu, dr_H, KR, KT);
		sprintf(filenameH2, "./CircularPlateData/H2_elastic-Nphi_%d-Npsi_%d-nu_%f-dr_%f-KR_%f-KT_%f.bin", NphiF, NpsiF, nu, dr_H, KR, KT);
	}
	
	while (!HLoaded){
		
		if (gs->getStormerVerletScheme()){
			ifstream fiH1(filenameH1, ios::in | ios::binary);
			if (fiH1.good()) {
					logMessage(1, "Loading H1 from cache");
					
					// found cache file, load it
					double *H1_aux = new double;
					int p;
					
						for (int i=0; i<Npsi; i++){
							p = 0;
							for (int j=0; j<Nphi; j++){
								for (int k=0; k<Nphi; k++){
									fiH1.seekg((i+j*NpsiF+k*NpsiF*NphiF)*sizeof(double), ios::beg);								
									fiH1.read((char *)H1_aux,sizeof(double));
									H1[i+p*Npsi] = *H1_aux * sqrt(e/2);
									p++; 	
								}
							}
						}
						HLoaded = true;
						fiH1.close();

						delete H1_aux;


			} 
			else {
				
				
				
				H_tensorCircular( NphiF, NpsiF, BC, nu, KR, KT, dr_H );
				
				logMessage(1, "Performed Htensor");
				
			}
		} 
		else {
			ifstream fiH0(filenameH0, ios::in | ios::binary);
			ifstream fiH1(filenameH1, ios::in | ios::binary);
			ifstream fiH2(filenameH2, ios::in | ios::binary);
			
			if (fiH0.good() && fiH1.good() && fiH2.good()) {
					logMessage(1, "Loading H tensors from cache");
					// found cache file, load it
					double *H0_aux = new double;
					double *H1_aux = new double;
					double *H2_aux = new double;
					int p;
					
						for (int i=0; i<Npsi; i++){
							p = 0;
							for (int j=0; j<Nphi; j++){
								for (int k=0; k<Nphi; k++){
									fiH0.seekg((i+j*NpsiF+k*NpsiF*NphiF)*sizeof(double), ios::beg);
									fiH0.read((char *)H0_aux,sizeof(double));
									
									fiH1.seekg((i+j*NpsiF+k*NpsiF*NphiF)*sizeof(double), ios::beg);
									fiH1.read((char *)H1_aux,sizeof(double));
																		
									fiH2.seekg((i+j*NpsiF+k*NpsiF*NphiF)*sizeof(double), ios::beg);
									fiH2.read((char *)H2_aux,sizeof(double));
									
									H0[i+p*Npsi] = *H0_aux * sqrt(e/2);
									H1[i+p*Npsi] = *H1_aux * sqrt(e/2);
									H2[i+p*Npsi] = *H2_aux * sqrt(e/2);
									

									p++; 	
								}
							}
						}
						HLoaded = true;

						fiH0.close();
						fiH1.close();
						fiH2.close();
						
						delete H0_aux;
						delete H1_aux;
						delete H2_aux;
			}
			else {
				H_tensorCircular( NphiF, NpsiF, BC, nu, KR, KT, dr_H );
				logMessage(1, "Performed Htensor");
				
			}
		}
	}
		
		
		
}




void CircularImperfectPlate::AxisymmetricCap(double H, double Rd, char ImperfectionType, double *Imperfection, int Nr, int Nth)
{
	
	int r, th;
	double dr;
	
	double y0;

	double R;
	
	double *y = new double [Nr];
	
	R = (H*H + Rd*Rd)/(2*H);
	
	dr = Rd/(Nr-1);
	

	switch (ImperfectionType){
		case 's': // spherical
			y0 = M_PI*((R*R*R-pow((R*R-Rd*Rd),1.5))*2/3-R*Rd*Rd);
			for (r = 0; r < Nr; r++){
				y[r]= -R + sqrt(R*R - r*dr*r*dr );
				y[r] = y[r]-y0/(M_PI*Rd*Rd);
			}
			break;
		case 'p': //Parabolic
			y0 = (-2*M_PI*H*pow(R, tau2 + 2))/(tau2 +2);
			for (r = 0; r < Nr; r++){
				y[r]= H*pow(r, tau2);
				y[r] = y[r]-y0/(M_PI*Rd*Rd);
			}
			break;
			
	}
	for (th = 0; th < Nth; th++){
		for (r = 0; r < Nr; r++){
			Imperfection[r*Nth+th] = y[r]; 
		}
	}
	ofstream fImp("./CircularPlateData/Imperfection.bin", ios::out | ios::binary);
	if (fImp.good()) {
	    fImp.write((const char *)Imperfection, Nr*Nth*sizeof(double));
	    fImp.close();
	}
	
	delete[] y;

}

void CircularImperfectPlate::ProjectionCoefficients(double *Imperfection, double *Ai, int Nphi, int modeType, double max_error,int Nr, int Nth, double nu, double KR, int *k, int *n, int *c, double *xkn){

	double zg;
	bool Include;
	double err_i, err_std, err_av;
	int i=0;
	int j;
	double norm=0;

	double *recons = new double[Nr*Nth];
	double *phi = new double[Nr*Nth];

	
	zg = 0; // 0 for axisymmetric caps

	
	for (j=0; j<(Nr*Nth); j++){
		recons[j]=zg;
		Imperfection[j] = Imperfection[j]/h;
	}

		
	err_i = 99999;
	
	
	while ((err_i>max_error)&&(i<Nphi)){

		Include = false;


		switch (modeType){
			case 0:
				Include = true;
				break;
			case 1: //Modes with k == 0	
				if (k[i]==0){
					Include = true;
				}
				break;
			case 2: // Modes with n == 0
				if (n[i]==0){
					Include = true;
				}
				break;
		}

		

		if (Include == true){

			
			ModeShape(phi, k[i], c[i], xkn[i], 1, nu, KR, Nr, Nth);
		
			norm = ScalarProductPolar(phi,phi,Nr,Nth);

			Ai[i] = ScalarProductPolar(Imperfection,phi,Nr,Nth)/sqrt(norm);

			err_std = 0;
			err_av = 0;
			for (j=0; j<(Nr*Nth); j++){
				recons[j] = recons[j] + Ai[i]*phi[j]/sqrt(norm); //phi should be normalized before computing Ai. Instead, it is divided by norm twice here.
				err_av = err_av + abs((Imperfection[j]-recons[j] + 1e-64)/(Imperfection[j] + 1e-64));		

			}
			
			err_av = err_av / (Nr*Nth);

			if ((Ai[i] != 0) || (err_av != 1 )){
			
				for (j=0; j<(Nr*Nth); j++){
					err_std = err_std + (abs((Imperfection[j]-recons[j] + 1e-64)/(Imperfection[j] + 1e-64)) - err_av)*(abs((Imperfection[j]-recons[j] + 1e-64)/(Imperfection[j] + 1e-64)) - err_av);
				}
				
				err_std = sqrt(err_std / (Nr*Nth - 1));

				if (err_i > err_av + err_std){
					err_i = err_av + err_std;
				}
			}

		} 
		else{
			Ai[i] = 0;
		}
		i++;
		
		

	}

	while (i<Nphi){
		Ai[i] = 0;
		i++;
	}

	delete[] recons;
	delete[] phi;

}

void CircularImperfectPlate::ModeShape(double *phi, int k,int c, double xkn, double Rd, double nu, double KR, int Nr, int Nth)
{

	double J0,J1,J2,I0,I1,I2,J,I;
	double Jtild,Itild, th_coef;
	double dx, dth;
	int x, th;
	double max=0;

	dx = xkn/(Nr-1);
	dth = 2*M_PI/(Nth-1);

    if (KR < 0){
        Jtild = gsl_sf_bessel_Jn(k,xkn);
        Itild = gsl_sf_bessel_In(k,xkn);
        
    	
    } 
    else {

		J2 = gsl_sf_bessel_Jn(k-2,xkn);
		J1 = gsl_sf_bessel_Jn(k-1,xkn);
		J0 = gsl_sf_bessel_Jn(k,xkn);
		I2 = gsl_sf_bessel_In(k-2,xkn);
		I1 = gsl_sf_bessel_In(k-1,xkn);
		I0 = gsl_sf_bessel_In(k,xkn);
	
		Jtild=xkn*xkn*J2+((nu-2*k+1)*xkn + KR)*J1 + (k*(k+1)*(1-nu) - KR*k)*J0;
		Itild=xkn*xkn*I2+((nu-2*k+1)*xkn + KR)*I1 + (k*(k+1)*(1-nu) - KR*k)*I0;
    }
    
	for (th = 0; th < Nth; th++){
		th_coef = cos(k*th*dth+(c-1)/2*M_PI);

		for (x = 0; x < Nr; x++){
			J = gsl_sf_bessel_Jn(k,x*dx);
			I = gsl_sf_bessel_In(k,x*dx);
			phi[x*Nth+th] = (J - (Jtild*I/Itild))*th_coef;
			
			if (max < phi[x*Nth+th]){
				max = phi[x*Nth+th];
			}
		}
	}

	for (x = 0; x<Nth*Nr; x++){
		phi[x]=phi[x]/max;
	}

}



double* CircularImperfectPlate::SortZeros( int xkLength, int Total, std::vector<CircularImperfectPlate::ModeXKN> xk, char BC)
{
	
	double * mode;
	int width = 6, j =0;
	
	mode = new double[Total*width];
	
	std::sort(xk.begin(), xk.begin()+xkLength);
	
	for (int i=0; i<xkLength; i++){
    	mode[j*width] = j;
    	mode[j*width+1] = xk[i].x;
    	mode[j*width+2] = (double)xk[i].k;
    	if ((BC=='f')&&(mode[j*width+2]!=1)){
    		mode[j*width+3] = (double)xk[i].n-1;
    	}else{
    		mode[j*width+3] = (double)xk[i].n;
    	}
    	
    	mode[j*width+4] = 1; 
    	mode[j*width+5] = mode[j*width+1]*mode[j*width+1];
    	if (mode[j*width+2]>0){
    		j++;
    		mode[j*width] = j;
			mode[j*width+1] = mode[(j-1)*width+1];
			mode[j*width+2] = mode[(j-1)*width+2];
			mode[j*width+3] = mode[(j-1)*width+3]; 
			mode[j*width+4] = 2; 
			mode[j*width+5] = mode[(j-1)*width+5];
    	}
    	
    	j++;
	}
	
	return mode;
}


double CircularImperfectPlate::CosCosCosIntegration(int k, int l, int m){
	double I = 0;
	
	if ((m == l+k)||(m == abs(l-k))){
		I = (M_PI/2);	
	}
	if (((k == 0) && (l == m)) || ((l == 0) && (k == m)) || ((m == 0) && (k == l))){
		
		I = M_PI;
	}
	if ((k == 0) && (l == 0) && (m == 0)){
		I = (2*M_PI);
	}
	return I;
}

double CircularImperfectPlate::CosSinSinIntegration(int k, int l, int m){
	double I = 0;
	
	if ((k == abs(l-m)) && (l != m) && (l != 0) && (m != 0)){
		I = (M_PI/2);
	}
	if ((k == (l+m)) && (l != 0) && (m != 0) ){
		I = (-M_PI/2);
	}
	if ((k == 0) && (l == m) && (l != 0)){
		I = M_PI;
	}
	return I;
}

double CircularImperfectPlate::HCoefficientCircular( int kp, int kq, int cp, int cq, double xip, double xiq, int ki, int ci, double zeta, double nu, double KR,  double dr_H ){
	double H = 0, rr;
	double JJ0, JJ1, dJJ0, ddJJ0, II0, II1, dII0, ddII0, Kkn, Lkn;
	double J, I, J2, J1, J0, I2, I1, I0, Jtild, Itild;

	int NH;
	NH = (int)(1/dr_H) + 1;
	
	double kk [2] = {kp, kq};
	double xi [2] = {xip, xiq};
	double k, xkn;
	
	double* R = new double[2*NH];
	double* dR = new double[2*NH];
	double* ddR = new double[2*NH];
	double* nRS = new double[NH];
	double* S = new double[NH];
	double* fctH1 = new double[NH];
	double* fctH2 = new double[NH];
	double h1, h2, beta1, beta2;
		
	
	if (!((ki == kp+kq) || (ki == abs(kp-kq)))){
		return H;
	}else if (((cp == cq) && (ci == 2)) || ((cp != cq) && (ci == 1))) {
		return H;
	}
	
	// Calculation of R, R', R'' (Transverse modes)
	
	for (int ii=0; ii<2; ii++){
		k = kk[ii];
		xkn = xi[ii];
		for(int nr=0; nr<NH; nr++){
			rr = nr*dr_H + 1e-32;
			
			JJ0 = gsl_sf_bessel_Jn(k,xkn*rr);
			JJ1 = gsl_sf_bessel_Jn(k+1,xkn*rr);
			dJJ0 = -JJ1*xkn + k*JJ0/rr;
			ddJJ0 = -(xkn*xkn + k/(rr*rr) - k*k/(rr*rr)) * JJ0 + xkn/rr*JJ1;
			
			II0 = gsl_sf_bessel_In(k,xkn*rr);
			II1 = gsl_sf_bessel_In(k+1,xkn*rr);
			dII0 = II1*xkn + k*II0/rr;
			ddII0 = (xkn*xkn - k/(rr*rr) + k*k/(rr*rr)) * II0 - xkn/rr*II1;
			
			
			
			if (KR < 0){
				J = gsl_sf_bessel_Jn(k,xkn);
				I = gsl_sf_bessel_In(k,xkn);
				
				R[ii*NH + nr] = I*JJ0 - J*II0;
				dR[ii*NH + nr] = I*dJJ0 - J*dII0;
				ddR[ii*NH + nr] = I*ddJJ0 - J*ddII0;
				
				nRS[nr] = R[ii*NH + nr]*R[ii*NH + nr];
				
			}else{
				J2 = gsl_sf_bessel_Jn(k-2,xkn);
				J1 = gsl_sf_bessel_Jn(k-1,xkn);
				J0 = gsl_sf_bessel_Jn(k,xkn);
				I2 = gsl_sf_bessel_In(k-2,xkn);
				I1 = gsl_sf_bessel_In(k-1,xkn);
				I0 = gsl_sf_bessel_In(k,xkn);
				
				Jtild=xkn*xkn*J2+((nu-2*k+1)*xkn + KR)*J1 + (k*(k+1)*(1-nu) - KR*k)*J0;
				Itild=xkn*xkn*I2+((nu-2*k+1)*xkn + KR)*I1 + (k*(k+1)*(1-nu) - KR*k)*I0;
				
				R[ii*NH + nr] = Itild*JJ0 - Jtild*II0;
				dR[ii*NH + nr] = Itild*dJJ0 - Jtild*dII0;
				ddR[ii*NH + nr] = Itild*ddJJ0 - Jtild*ddII0;
				nRS[nr] = R[ii*NH + nr]*R[ii*NH + nr];
			}
			
			
		}
		
		
			
		Kkn = sqrt(1/Trapz1DPolar( nRS, 1e-32, 1, NH));
		
		
		
		if (k == 0){
			Kkn = Kkn/sqrt(2*M_PI);
		}else{
			Kkn = Kkn/sqrt(M_PI);
		}
		
		
		
		for(int nr=0; nr<NH; nr++){
			R[ii*NH + nr] = Kkn*R[ii*NH + nr];
			dR[ii*NH + nr] = Kkn*dR[ii*NH + nr];
			ddR[ii*NH + nr] = Kkn*ddR[ii*NH + nr];
			
		}
		
	}
	
	// Calculation of S, S', S'' (In-plane modes)
	
	if (KR<0){
		J2 = gsl_sf_bessel_Jn(ki-2,zeta);
		J1 = gsl_sf_bessel_Jn(ki-1,zeta);
		J0 = gsl_sf_bessel_Jn(ki,zeta);
		I2 = gsl_sf_bessel_In(ki-2,zeta);
		I1 = gsl_sf_bessel_In(ki-1,zeta);
		I0 = gsl_sf_bessel_In(ki,zeta);
		
		Jtild = zeta*zeta*J2 + (-nu-2*ki+1)*zeta*J1 + (ki*(ki+1) + nu*ki*(1-ki))*J0;
		Itild = zeta*zeta*I2 + (-nu-2*ki+1)*zeta*I1 + (ki*(ki+1) + nu*ki*(1-ki))*I0;
		
		for(int nr=0; nr<NH; nr++){
			rr = nr*dr_H + 1e-32;
			J = gsl_sf_bessel_Jn(ki,zeta*rr);
			I = gsl_sf_bessel_In(ki,zeta*rr); 
			
			S[nr] = Itild*J - Jtild*I;
			
			nRS[nr] = S[nr]*S[nr];
			
			
		}
	} 
	else{
		J = gsl_sf_bessel_Jn(ki,zeta);
		I = gsl_sf_bessel_In(ki,zeta); 
		for(int nr=0; nr<NH; nr++){
			rr = nr*dr_H + 1e-32;
			JJ0 = gsl_sf_bessel_Jn(ki,zeta*rr);
			II0 = gsl_sf_bessel_In(ki,zeta*rr); 
			
			S[nr] = I*JJ0 - J*II0;
			
			
			nRS[nr] = S[nr]*S[nr];
								
			
		}
	}
	
	Lkn = sqrt(1/Trapz1DPolar( nRS, 1e-32, 1, NH));
							
	if (ki == 0){
		Lkn = Lkn/sqrt(2*M_PI);
	}else{
		Lkn = Lkn/sqrt(M_PI);
	}
	
	for(int nr=0; nr<NH; nr++){
		S[nr] = Lkn*S[nr];
	}
	
	// Computation of the coefficients
	
	for(int nr=0; nr<NH; nr++){
		rr = nr*dr_H + 1e-32;
		//fctH1 = S.*(ddR(1,:).*(dR(2,:)-kq^2*R(2,:)./rr)+ddR(2,:).*(dR(1,:)-kp^2*R(1,:)./rr));
		fctH1[nr] = S[nr]*(ddR[nr]*(dR[1*NH + nr] - kq*kq*R[1*NH + nr]/rr) + ddR[1*NH + nr]*(dR[nr] - kp*kp*R[nr]/rr));
		//fctH2 = S.*(dR(1,:)-R(1,:)./rr).*(dR(2,:)-R(2,:)./rr)./rr;
		fctH2[nr] = S[nr]*(dR[nr] - R[nr]/rr)*(dR[1*NH+nr]-R[1*NH+nr]/rr)/rr;	
	}
	
	// Radius dependent term
	h1 = Trapz1DCartesian( fctH1, 1e-32, 1, NH);
	h2 = Trapz1DCartesian( fctH2, 1e-32, 1, NH);
	
	// Theta dependent term
	if (cp ==1){
		if (cq == 1){
			beta1 = CosCosCosIntegration(kp,kq,ki);
			beta2 = CosSinSinIntegration(ki,kp,kq);
		} else{
			 beta1 = CosSinSinIntegration(kp,kq,ki);
			 beta2 = -CosSinSinIntegration(kq,kp,ki);
		}
	} else {
		if (cq==1){
			beta1 = CosSinSinIntegration(kq,kp,ki);
			beta2 = -CosSinSinIntegration(kp,kq,ki);
		} else{
			beta1 = CosSinSinIntegration(ki,kp,kq);
			beta2 = CosCosCosIntegration(kq,kp,ki);
		}	
	}
	
	H = h1*beta1 -2*kp*kq*h2*beta2;
	
	return H;
}

void CircularImperfectPlate::H_tensorCircular(int NphiF, int NpsiF, char BC, double nu, double KR, double KT, double dr_H )
{
	
	char filenameH0[1000], filenameH1[1000], filenameH2[1000];//, filenameT[1000], filenameL[1000];
	
	if (BC == 'f'){
			sprintf(filenameH0, "./CircularPlateData/H0_free-Nphi_%d-Npsi_%d-nu_%f-dr_%f.bin", NphiF, NpsiF, nu, dr_H);
			sprintf(filenameH1, "./CircularPlateData/H1_free-Nphi_%d-Npsi_%d-nu_%f-dr_%f.bin", NphiF, NpsiF, nu, dr_H);
			sprintf(filenameH2, "./CircularPlateData/H2_free-Nphi_%d-Npsi_%d-nu_%f-dr_%f.bin", NphiF, NpsiF, nu, dr_H);
			//sprintf(filenameT, "./CircularPlateData/mode_t_free-nu_%f.bin",nu);
			//sprintf(filenameL, "./CircularPlateData/mode_l_free.bin");
			KR = 0;
	} else if (BC == 'c'){
			sprintf(filenameH0, "./CircularPlateData/H0_clamped-Nphi_%d-Npsi_%d-nu_%f-dr_%f.bin", NphiF, NpsiF, nu, dr_H);
			sprintf(filenameH1, "./CircularPlateData/H1_clamped-Nphi_%d-Npsi_%d-nu_%f-dr_%f.bin", NphiF, NpsiF, nu, dr_H);
			sprintf(filenameH2, "./CircularPlateData/H2_clamped-Nphi_%d-Npsi_%d-nu_%f-dr_%f.bin", NphiF, NpsiF, nu, dr_H);
			//sprintf(filenameT, "./CircularPlateData/mode_t_clamped.bin");
			//sprintf(filenameL, "./CircularPlateData/mode_l_clamped-nu_%f.bin", nu);
			KR = -1;
	} else {
			sprintf(filenameH0, "./CircularPlateData/H0_elastic-Nphi_%d-Npsi_%d-nu_%f-dr_%f-KR_%f-KT_%f.bin", NphiF, NpsiF, nu, dr_H, KR, KT);
			sprintf(filenameH1, "./CircularPlateData/H1_elastic-Nphi_%d-Npsi_%d-nu_%f-dr_%f-KR_%f-KT_%f.bin", NphiF, NpsiF, nu, dr_H, KR, KT);
			sprintf(filenameH2, "./CircularPlateData/H2_elastic-Nphi_%d-Npsi_%d-nu_%f-dr_%f-KR_%f-KT_%f.bin", NphiF, NpsiF, nu, dr_H, KR, KT);

	}
	
	double dx = 1e-3, xmax = 100;
	int *ktF, *ctF, *klF, *clF;
	double *xiF, *ovF, *zlF;
	
    xiF = new double[NphiF];
	ktF = new int[NphiF];
	ctF = new int[NphiF];
	ovF = new double[NphiF];
	zlF = new double[NpsiF];
	klF = new int[NpsiF];
	clF = new int[NpsiF];
	
    TransverseModeCalculator(dx, xmax, BC, nu, KR, KT, NphiF, xiF, ktF, ctF, ovF);    
    InPlaneModeCalculator(dx, xmax, BC, nu, NpsiF, zlF, klF, clF);
    
	double *H0_aux = new double[NpsiF*NphiF*NphiF];
	double *H1_aux = new double[NpsiF*NphiF*NphiF];
	double *H2_aux = new double[NpsiF*NphiF*NphiF];

    for (int p=0; p<NphiF; p++){
    	for (int q=p; q<NphiF; q++){
    		for (int i = 0; i<NpsiF; i++){
    			
                H0_aux[i + p*NpsiF + q*NpsiF*NphiF] = HCoefficientCircular( ktF[p], ktF[q], ctF[p], ctF[q], xiF[p], xiF[q], klF[i], clF[i], zlF[i], nu, KR, dr_H );
                
                H0_aux[i + q*NpsiF + p*NpsiF*NphiF] = H0_aux[i + p*NpsiF + q*NpsiF*NphiF];
                
                H1_aux[i + p*NpsiF + q*NpsiF*NphiF] = H0_aux[i + p*NpsiF + q*NpsiF*NphiF]/zlF[i] / zlF[i];
                
                H1_aux[i + q*NpsiF + p*NpsiF*NphiF] = H0_aux[i + q*NpsiF + p*NpsiF*NphiF]/zlF[i] / zlF[i];
                
                H2_aux[i + p*NpsiF + q*NpsiF*NphiF] = H0_aux[i + p*NpsiF + q*NpsiF*NphiF]/zlF[i] / zlF[i] / zlF[i] / zlF[i];
                
                H2_aux[i + q*NpsiF + p*NpsiF*NphiF] = H0_aux[i + q*NpsiF + p*NpsiF*NphiF]/zlF[i] / zlF[i] / zlF[i] / zlF[i];
                
            
    		}
    	}
    }
    
	ofstream foH0(filenameH0, ios::out | ios::binary);
	if (foH0.good()) {
	    foH0.write((const char *)H0_aux, NphiF*NphiF*NpsiF*sizeof(double));
	    foH0.close();
	}
	
	ofstream foH1(filenameH1, ios::out | ios::binary);
	if (foH1.good()) {
	    foH1.write((const char *)H1_aux, NphiF*NphiF*NpsiF*sizeof(double));
	    foH1.close();
	}
	
	ofstream foH2(filenameH2, ios::out | ios::binary);
	if (foH2.good()) {
	    foH2.write((const char *)H2_aux, NphiF*NphiF*NpsiF*sizeof(double));
	    foH2.close();
	}
	
	delete[] H0_aux;
	delete[] H1_aux;
	delete[] H2_aux;
	

}


double CircularImperfectPlate::Norm_Modes(int k, double xkn, double R){
	
	int i = 0;
	double JJ0, II0, J2, J1, J0, I2, I1, I0, Jkn, Ikn, Jtild, Itild;
	
	double *Rkn = new double[Nr];
	double Kkn = 0;
	
	double dr = R/(Nr-1);

	if ((BC == 'f') || (BC == 'e')){
		
		J2 = gsl_sf_bessel_Jn(k-2,xkn);
		J1 = gsl_sf_bessel_Jn(k-1,xkn);
		J0 = gsl_sf_bessel_Jn(k,xkn);
		I2 = gsl_sf_bessel_In(k-2,xkn);
		I1 = gsl_sf_bessel_In(k-1,xkn);
		I0 = gsl_sf_bessel_In(k,xkn);
		
		Jtild = xkn*xkn*J2+((nu-2*k+1)*xkn + KR)*J1 + (k*(k+1)*(1-nu) - KR*k)*J0;
		Itild = xkn*xkn*I2+((nu-2*k+1)*xkn + KR)*I1 + (k*(k+1)*(1-nu) - KR*k)*I0;

		
		for (i=0; i<Nr; i++){

			JJ0 = gsl_sf_bessel_Jn(k,xkn*i*dr);
			II0 = gsl_sf_bessel_In(k,xkn*i*dr);

			Rkn[i] = JJ0 - ((Jtild*II0)/Itild);

			Rkn[i] = Rkn[i]*Rkn[i];
			
			

			
		}
	}
	else {
		Jkn = gsl_sf_bessel_Jn(k,xkn);
		Ikn = gsl_sf_bessel_In(k,xkn);
		
		for (i=0; i<Nr; i++){

			JJ0 = gsl_sf_bessel_Jn(k,xkn*i*dr);
			II0 = gsl_sf_bessel_In(k,xkn*i*dr);

			Rkn[i] = (JJ0*Ikn - Jkn*II0);
			Rkn[i] = Rkn[i]*Rkn[i];

			
		}
	}
		
	Kkn = Trapz1DPolar(Rkn, 0, 1, Nr);
		
	Kkn = 1 / (sqrt(Kkn*M_PI));

	if (k == 0){
		Kkn = Kkn/sqrt(2);
	}
	
	delete[] Rkn;

	return Kkn;
	
}




