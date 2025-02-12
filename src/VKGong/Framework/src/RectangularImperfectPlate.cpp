/*
 * NeSS Framework Code
 *
 * Copyright (c) The University of Edinburgh, 2016. All rights reserved.
 *
 * Author: James Perry (j.perry@epcc.ed.ac.uk)
 */

#include "RectangularImperfectPlate.h"
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
using namespace std;

void getCacheDirectory(char *buf, int bufsize);

RectangularImperfectPlate::RectangularImperfectPlate(string name, int iNphi, int iNpsi, double Lx, double Ly, double h, double H, char ImperfectionType,
		double xWidth, double yWidth, int modeType, double nu, double Young, double rho, char BC, int Nx, int Ny, double dFac, double dExp, double dCons, int fs)
    : Component(name)
{
    int i;
    GlobalSettings *gs = GlobalSettings::getInstance();

    logMessage(1, "Creating RectangularImperfectPlate: modes=%d,%d, size=%fx%f, h=%f, nu=%f, E=%f, rho=%f", iNphi, iNpsi, Lx, Ly, h, nu, Young, rho);
    logMessage(1, "Imperfection: H = %f, ImperfectionType = %c, xWidth = %f, yWidth = %f", H, ImperfectionType, xWidth, yWidth);
    

    this->Lx = Lx;
    this->Ly = Ly;
    
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
    
    this->ImperfectionType = ImperfectionType;
    this->modeType = modeType;
    this->BC = BC;

    this->Nx = Nx;
	this->Ny = Ny;
	
	this->xWidth = xWidth;
	this->yWidth = yWidth;
    
	this->e = (Lx*Ly/4)/rho*Young;

    // allocate state arrays
    u = new double[Nphi];
    u1 = new double[Nphi];
    u2 = new double[Nphi];

    memset(u, 0, Nphi * sizeof(double));
    memset(u1, 0, Nphi * sizeof(double));
    memset(u2, 0, Nphi * sizeof(double));
    
    int S;
    double *coeff0, *coeff1, *coeff2;

    ov = new double [Nphi];
    kx = new int [Nphi];
    ky = new int [Nphi];
    
    
    
    // allocate temporary buffers corresponding to every scheme.
    if (gs->getStormerVerletScheme()){
    	// perform eigencalc
   	
    	  	
		coeff1 = new double[Npsi*Npsi*Npsi*Npsi*sizeof(double)];

		coeff0 = NULL;
		coeff2 = NULL;
		
		S = AiryStressFactorsCalculation(Npsi, Lx, Ly, coeff0, coeff1, coeff2);
		

								
		logMessage(1, "Performed AiryStressFactorsCalculation, S=%d", S);
		
		
		t0 = new double[S*Nphi];
		t1 = new double[S*Nphi];
		t2 = new double[S];
		t3 = new double[S];
		qAi = new double[Nphi];
		q2Ai = new double[Nphi];
		G = new double[Nphi];

		
		memset(G, 0, Nphi * sizeof(double));
		
		
		
		H0 = NULL;
		H2 = NULL;

		H1 = new double[Nphi*Nphi*S];



    }
    else{
    	

    	
		// perform eigencalc
		coeff0 = new double[Npsi*Npsi*Npsi*Npsi*sizeof(double)];
		coeff1 = new double[Npsi*Npsi*Npsi*Npsi*sizeof(double)];
		coeff2 = new double[Npsi*Npsi*Npsi*Npsi*sizeof(double)];
		
		S = AiryStressFactorsCalculation(Npsi, Lx, Ly, coeff0, coeff1, coeff2);
				    			    		
		logMessage(1, "Performed AiryStressFactorsCalculation, S=%d", S);
		
    	t0 = new double[S*Nphi];
    	t0t = new double[S*Nphi];
		t1 = new double[S*Nphi];
		t4 = new double[S];
		t5 = new double[S*Nphi];
		t7 = new double[S*Nphi];
		t8 = new double[S*Nphi];
		qAi = new double[Nphi];
		qq1 = new double[Nphi];
		
		eta = new double[S];
		eta1 = new double[S];
		eta2 = new double[S];
		eta_temp = new double[S];
		
		memset(eta, 0, S * sizeof(double));
		memset(eta1, 0, S * sizeof(double));
		memset(eta2, 0, S * sizeof(double));
		memset(eta_temp, 0, S * sizeof(double));
		
		Ga = new double[Nphi];
		Gb = new double[Nphi];
		G = new double[Nphi];
		
		memset(G, 0, Nphi * sizeof(double));
		memset(Ga, 0, Nphi * sizeof(double));
		memset(Gb, 0, Nphi * sizeof(double));
		
		Gc1 = new double[S];
		Gc2 = new double[S];
		Gc = new double[S];
		
		memset(Gc, 0, S * sizeof(double));
		memset(Gc1, 0, S * sizeof(double));
		memset(Gc2, 0, S * sizeof(double));
		
		G0 = new double[Nphi*Nphi];
		mat_imp = new double[Nphi*Nphi];
		
		memset(G0, 0, Nphi * Nphi * sizeof(double));
		memset(mat_imp, 0, Nphi * Nphi * sizeof(double));

		
		H0 = new double[Nphi*Nphi*S];
		H1 = new double[Nphi*Nphi*S];
		H2 = new double[Nphi*Nphi*S];
		
		
    }
    
    ComputeTransverseEigenfrequenciesRectangular(Nphi, Lx, Ly, h, Young, rho, nu, ov, kx, ky);
      
    
    LoadHTensor(coeff0, coeff1, coeff2, Nphi, Npsi, Lx, Ly, S, H0, H1, H2);
    

    //this->Npsi = S;
	//Npsi = S;

	if (coeff0) {delete[] coeff0;}
	if (coeff2) {delete[] coeff2;}
	delete[] coeff1;
    
    Ai = new double[Nphi];
	memset(Ai, 0, Nphi * sizeof(double));
	
	if (H>0){

	// see if cached Ai for this imperfection exists
	char cachefile[1000];
	sprintf(cachefile, "./RectangularPlateData/Imperfection_Ai-Nphi_%d-H_%f-Lx_%f-Ly_%f_%c.bin", Nphi, H, Lx, Ly, ImperfectionType); //[angels] Aclarir el format per guardar el tensor H
	ifstream fi2(cachefile, ios::in | ios::binary);
	if (fi2.good()) {
		logMessage(1, "Loading Ai from file.");

	// found cache file, load it
	   fi2.read((char *)Ai, Nphi*sizeof(double));
	   fi2.close();
    }
	else {
	   double *Imperfection = new double[Nx*Ny*sizeof(double)];  


	   RectangularImperfection( Imperfection, Lx, Ly, H, Nx, Ny, ImperfectionType, xWidth, yWidth);

	   
	   double max_error = 0.1;
   
	   ProjectionCoefficients(Imperfection, Ai, Nphi, modeType, max_error, Nx, Ny, Lx, Ly);
	   
			   
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

    memset(Cnorm, 1, Nphi * sizeof(double));
    
    printf("plate fs = %d", fs);
    
    plate_def(&fs);

    logMessage(1, "Performed plate_def, sample rate=%d", fs);


    gs->setSampleRate((double)fs);



    logMessage(1, "Final mode counts: Npsi=%d, Nphi=%d", this->Npsi, this->Nphi);


    




}

RectangularImperfectPlate::~RectangularImperfectPlate()
{

    delete[] Cnorm;
    delete[] ov;
    delete[] kx;
    delete[] ky;

}

void RectangularImperfectPlate::plate_def(int *fs)
{
    double fs_lim;
    double fsd;
    double xi;
    int m;
    double *c = new double[Nphi];

    GlobalSettings *gs = GlobalSettings::getInstance();
    
    /* calculate sampling rate */
    fs_lim = ov[(Nphi-1)] / 2.0;
    
    if (*fs < fs_lim) {    
		*fs = (int)(fs_lim * 3); // fs is stored in dimensional form
		printf("Warning: The sampling rate introduced by user is too low. It will be modified to fulfill the stability limit. \n");
	}
    
    printf("Simulation fs = %d \n", *fs);
    
    fsd = (double)*fs;

    /* damping ratios */
    for (m = 0; m < Nphi; m++) {
    	c[m] = dFac * pow(ov[m], dExp) + dCons;    	    	
    }

    for (m = 0; m < Nphi; m++) {
	C[m] = (fsd*fsd + c[m]*fsd);
	C1[m] = (-2.0*fsd*fsd + ov[m]*ov[m]);
	C2[m] = (fsd*fsd - c[m]*fsd);
	
	
	if (gs->getStormerVerletScheme()){
		Cnorm[m] = C[m]; /* Cnorm is just a copy of C  used if StormerVerlet*/
		C1[m] = C1[m] / C[m];
		C2[m] = C2[m] / C[m];
	}
	else {
		Cnorm[m] = 1;
	}
    }
    delete[] c;
}



#define PI2 (M_PI*M_PI)
#define PI3 (M_PI*M_PI*M_PI)
#define PI4 (M_PI*M_PI*M_PI*M_PI)
#define PI5 (M_PI*M_PI*M_PI*M_PI*M_PI)

void RectangularImperfectPlate::H_tensorRectangular(double *coeff0, double *coeff1, double *coeff2, int Nphi, int Npsi, double Lx,
				double Ly, int S, double *H0, double *H1, double *H2)
{
    int m, n, p;
    int n1, p1, n2, p2;
    double *mthing;

    double tmp, tmp0, tmp1, tmp2;
    double m1, m2, m3, m4, m5, m6;
    double fac, rowmax0, rowmax1, rowmax2;

    //double *mode_t;
    
    char filenameH0[1000], filenameH1[1000], filenameH2[1000];

    /* allocate temporary storage */
    mthing = new double[S*Nphi*Nphi];
   // mode_t = new double[Nphi*Nphi*3];

    /* generate mode_t */
    //ComputeTransverseEigenfrequenciesRectangular(Nphi, Lx, Ly, ov, kx, ky );


    
    /* generate mthing (column major) */
    for (m = 0; m < S; m++) {
    
	int ma = m / Npsi; /* needed for g1, g2, g5 */
	int mb = m % Npsi; /* needed for g3, g4, g6 */

    for (n = 0; n < Nphi; n++) {
	    n1 = kx[n];
	    n2 = ky[n];
	    for (p = 0; p < Nphi; p++) {
		p1 = kx[p];
		p2 = ky[p];

		

		tmp = (i1_mat(Npsi,Nphi,Lx,ma+1,n1,p1) + i2_mat(Npsi,Nphi,Lx,ma+1,n1,p1) +
		       i3_mat(Npsi,Nphi,Lx,ma+1,n1,p1) + i4_mat(Npsi,Nphi,Lx,ma+1,n1,p1) +
		       i5_mat(Npsi,Nphi,Lx,ma+1,n1,p1));
		m1 = tmp * ((double)(n1*n1));
		m2 = tmp * ((double)(p1*p1));
		
		m5 = (i9_mat (Npsi,Nphi,Lx,ma+1,n1,p1) + i10_mat(Npsi,Nphi,Lx,ma+1,n1,p1) +
		      i11_mat(Npsi,Nphi,Lx,ma+1,n1,p1) + i12_mat(Npsi,Nphi,Lx,ma+1,n1,p1) +
		      i13_mat(Npsi,Nphi,Lx,ma+1,n1,p1)) * (double)p1 * (double)n1;

		tmp = (i1_mat(Npsi,Nphi,Ly,mb+1,n2,p2) + i2_mat(Npsi,Nphi,Ly,mb+1,n2,p2) +
		       i3_mat(Npsi,Nphi,Ly,mb+1,n2,p2) + i4_mat(Npsi,Nphi,Ly,mb+1,n2,p2) +
		       i5_mat(Npsi,Nphi,Ly,mb+1,n2,p2));
		m3 = tmp * ((double)(n2*n2));
		m4 = tmp * ((double)(p2*p2));

		m6 = (i9_mat (Npsi,Nphi,Ly,mb+1,n2,p2) + i10_mat(Npsi,Nphi,Ly,mb+1,n2,p2) +
		      i11_mat(Npsi,Nphi,Ly,mb+1,n2,p2) + i12_mat(Npsi,Nphi,Ly,mb+1,n2,p2) +
		      i13_mat(Npsi,Nphi,Ly,mb+1,n2,p2)) * (double)p2 * (double)n2;
		

		mthing[m + (n*S) + (p*Nphi*S)] = m1*m4 + m2*m3 - 2.0*m5*m6;
		
		

		
	    }
	}
    }
    
    

    fac = 4.0 * PI4 / (Lx*Lx*Lx) / (Ly*Ly*Ly);

    /* generate next row of H1 */
    for (n = 0; n < S; n++) {
	rowmax0 = 0.0;
	rowmax1 = 0.0;
	rowmax2 = 0.0;

	/*
	 * this row is product of coeff1 column n, with:
	 *   (m1.*m4 + m2.*m3 - 2*m5.*m6): mthing
	 */
	for (p = 0; p < (Nphi*Nphi); p++) {
	    tmp0 = 0.0;
	    tmp1 = 0.0;
	    tmp2 = 0.0;
	    for (m = 0; m < S; m++) {
		tmp1 += coeff1[(n*S)+m] * mthing[(p*S)+m];
		if (H0 != NULL){
			tmp0 += coeff0[(n*S)+m] * mthing[(p*S)+m];
			tmp2 += coeff2[(n*S)+m] * mthing[(p*S)+m];
		}
		
		
	    }

	    /* apply a scalar factor */
	    tmp1 *= fac;
	    if (H0 != NULL){
	    	tmp0 *= fac;
	        tmp2 *= fac;
	    }

	    /* store in H1, in column-major order */
	    H1[(p*S)+n] = tmp1;
	        
	    
	    if (H0 != NULL){
	    	H0[(p*S)+n] = tmp0;
	    	H2[(p*S)+n] = tmp2;
	    }
	    
	   

	    /* find the row max */
	    if (fabs(tmp1) > rowmax1) rowmax1 = fabs(tmp1);
	    if (H0 != NULL){
	    	if (fabs(tmp0) > rowmax0) rowmax0 = fabs(tmp0);
	    	if (fabs(tmp2) > rowmax2) rowmax2 = fabs(tmp2);
	    }
	}

	/* now zero all values that are 10 orders of magnitude lower than the row max */
	for (p = 0; p < (Nphi*Nphi); p++) {
	    if ((fabs(H1[(p*S)+n]) / rowmax1) < 1e-10) H1[(p*S)+n] = 0.0;
		    
	    if (H0 != NULL){
	    	if ((fabs(H0[(p*S)+n]) / rowmax0) < 1e-10) H0[(p*S)+n] = 0.0;
	    	if ((fabs(H2[(p*S)+n]) / rowmax2) < 1e-10) H2[(p*S)+n] = 0.0;
	    }
	}
    }
    
    
	sprintf(filenameH1, "./RectangularPlateData/H1_SS-Nphi_%d-S_%d-Lx_%f-Ly_%f.bin", Nphi, S, Lx, Ly);
	ofstream foH1(filenameH1, ios::out | ios::binary);
	if (foH1.good()) {
	    foH1.write((const char *)H1, Nphi*Nphi*S*sizeof(double));
	    foH1.close();
	}
	
	if (H0 != NULL){
		sprintf(filenameH0, "./RectangularPlateData/H0_SS-Nphi_%d-S_%d-Lx_%f-Ly_%f.bin", Nphi, S, Lx, Ly);	
		sprintf(filenameH2, "./RectangularPlateData/H2_SS-Nphi_%d-S_%d-Lx_%f-Ly_%f.bin", Nphi, S, Lx, Ly);
		
		ofstream foH0(filenameH0, ios::out | ios::binary);
		if (foH0.good()) {
			foH0.write((const char *)H0, Nphi*Nphi*S*sizeof(double));
			foH0.close();
		}

		ofstream foH2(filenameH2, ios::out | ios::binary);
		if (foH2.good()) {
			foH2.write((const char *)H2, Nphi*Nphi*S*sizeof(double));
			foH2.close();
		}
	}
    
    
    

    delete[] mthing;
    //delete[] mode_t;
}


static int mode_t_compare(const void *i1, const void *i2)
{
    double *mode_t1 = (double *)i1;
    double *mode_t2 = (double *)i2;
    if (mode_t1[0] < mode_t2[0]) return -1;
    if (mode_t1[0] > mode_t2[0]) return 1;

    if (mode_t1[1] < mode_t2[1]) return -1;
    if (mode_t1[1] > mode_t2[1]) return 1;

    if (mode_t1[2] < mode_t2[2]) return -1;
    if (mode_t1[2] > mode_t2[2]) return 1;

    return 0;
}

/* y should be big enough to hold a (Nphi*Nphi)x3 array, row major */
void RectangularImperfectPlate::ComputeTransverseEigenfrequenciesRectangular(int Nphi, double Lx, double Ly, double h, double Young, double rho, double nu, double* ov, int* kx, int* ky)
{
    int m, n, ii;
    double *mode_t;
    
    mode_t = new double[Nphi*Nphi*3];
    
    double D = Young * (h*h*h) / 12.0 / (1.0 - (nu*nu));

    ii = 0;
    for (m = 1; m <= Nphi; m++) {
	for (n = 1; n <= Nphi; n++) {
	    mode_t[(ii*3)+0] = sqrt(D/rho/h) * ((((double)m)*M_PI/Lx)*(((double)m)*M_PI/Lx) +
		(((double)n)*M_PI/Ly)*(((double)n)*M_PI/Ly));
	    mode_t[(ii*3)+1] = (double)m;
	    mode_t[(ii*3)+2] = (double)n;
	    
	    ii++;
	}
    }
    qsort(mode_t, Nphi*Nphi, 3*sizeof(double), mode_t_compare);
    
    for (ii = 0; ii<Nphi; ii++){
    	ov[ii] = mode_t[(ii*3)+0];
    	kx[ii] = mode_t[(ii*3)+1];
    	ky[ii] = mode_t[(ii*3)+2];
    	
    }
    delete[] mode_t;
    
}

/* m, n, p are one-based Matlab indices */
double RectangularImperfectPlate::i1_mat(int Npsi, int Nphi, double L, int m, int n, int p)
{
    int m1 = m - 1;
    if ((m1 == 0) && (n == p)) {
	return L/2.0;
    }
    else if ((m1 == (p-n)) || (m1 == (n-p))) {
	return L/4.0;
    }
    else if ((m1 == (-n-p)) || (m1 == (n+p))) {
	return -L/4.0;
    }
    return 0.0;
}

double RectangularImperfectPlate::i2_mat(int Npsi, int Nphi, double L, int m, int n, int p)
{
    int m1 = m - 1;
    double L4 = L*L*L*L;
    double L5 = L4*L;
    double dp = (double)p;
    double p3 = p*p*p;
    double p5 = p3*p*p;
    double npp = (double)n + (double)p;
    double npp4 = npp*npp*npp*npp;
    double npp5 = npp4*npp;
    double nmp = (double)n - (double)p;
    double nmp4 = nmp*nmp*nmp*nmp;
    double nmp5 = nmp4*nmp;
    double m1pm1 = 1.0;
    if (m1 & 1) m1pm1 = -1.0;

    if (n == p) {
	return (+15.0/L4*(m1pm1 + 1.0)) * (L5*(4.0*PI5*p5 - 20.0*PI3*p3 + 30.0*M_PI*dp)) /
	    (40.0*PI5*p5);
    }
    return -(+15.0/L4*(m1pm1 + 1.0)) * (8796093022208.0*L*((sin(M_PI*npp)*((1713638851887625.0*L4*npp4)/17592186044416.0 - (8334140006820045.0*L4*npp*npp)/70368744177664.0 + 24.0*L4) + 4.0*M_PI*L*L*cos(M_PI*npp)*npp*((2778046668940015.0*L*L*npp*npp)/281474976710656.0 - 6.0*L*L))/npp5 - (sin(M_PI*nmp)*((1713638851887625.0*L4*nmp4)/17592186044416.0 - (8334140006820045.0*L4*nmp*nmp)/70368744177664.0 + 24.0*L4) + 4.0*M_PI*L*L*cos(M_PI*nmp)*nmp*((2778046668940015.0*L*L*nmp*nmp)/281474976710656.0 - 6.0*L*L))/nmp5))/5383555227996211.0;
}

double RectangularImperfectPlate::i3_mat(int Npsi, int Nphi, double L, int m, int n, int p)
{
    int m1 = m - 1;
    double L3 = L*L*L;
    double L4 = L3*L;
    double dp = (double)p;
    double p4 = dp*dp*dp*dp;
    double nmp = (double)n - (double)p;
    double nmp4 = nmp*nmp*nmp*nmp;
    double npp = (double)n + (double)p;
    double npp4 = npp*npp*npp*npp;
    double m1pm1 = 1.0;
    if (m1 & 1) m1pm1 = -1.0;

    if (n == p) {
	return -(-4.0/L3*(7.0*m1pm1 + 8.0))*(L4*(6.0*PI2*dp*dp - 2.0*PI4*p4))/(16.0*PI4*p4);
    }
    return (-4.0/L3*(7.0*m1pm1 + 8.0))*(L*((6.0*L3)/nmp4 - (6.0*L3)/npp4))/(2.0*PI4) + (-4.0/L3*(7.0*m1pm1 + 8.0))*(L*(3.0*L*cos(M_PI*npp)*(2.0*L*L - L*L*PI2*npp*npp)/npp4 - (3.0*L*cos(M_PI*nmp)*(2.0*L*L - L*L*PI2*nmp*nmp))/nmp4))/(2.0*PI4);
}

double RectangularImperfectPlate::i4_mat(int Npsi, int Nphi, double L, int m, int n, int p)
{
    int m1 = m - 1;
    double L3 = L*L*L;
    double dp = (double)p;
    double p3 = dp*dp*dp;
    double nmp = (double)n - (double)p;
    double npp = (double)n + (double)p;
    double m1pm1 = 1.0;
    if (m1 & 1) m1pm1 = -1.0;

    if (n == p) {
	return -(6.0/(L*L)*(2.0*m1pm1 + 3.0)) * (L3*(6.0*M_PI*dp - 4.0*PI3*p3)) / (24.0*PI3*p3);
    }
    return (6/(L*L)*(2.0*m1pm1 + 3.0)) * (L3*cos(M_PI*nmp)) / (PI2*nmp*nmp) - (6.0/(L*L)*(2.0*m1pm1 + 3.0)) * (L3*cos(M_PI*npp)) / (PI2*npp*npp);
}

double RectangularImperfectPlate::i5_mat(int Npsi, int Nphi, double L, int m, int n, int p)
{
    if (n == p) return -L/2.0;
    return 0.0;
}

double RectangularImperfectPlate::i9_mat(int Npsi, int Nphi, double L, int m, int n, int p)
{
    int m1 = m - 1;
    if ((m1 == 0) && (n == p)) return L/2.0;
    else if ((m1 == (p-n)) || (m1 == (n-p)) || (m1 == (-n-p)) || (m1 == (n+p))) return L/4.0;
    return 0.0;
}

double RectangularImperfectPlate::i10_mat(int Npsi, int Nphi, double L, int m, int n, int p)
{
    int m1 = m - 1;
    double L4 = L*L*L*L;
    double L5 = L4*L;
    double dn = (double)n;
    double n3 = dn*dn*dn;
    double n5 = n3*dn*dn;
    double npp = (double)n + (double)p;
    double npp2 = npp*npp;
    double npp4 = npp2*npp2;
    double nmp = (double)n - (double)p;
    double nmp2 = nmp*nmp;
    double nmp4 = nmp2*nmp2;
    double m1pm1 = 1.0;
    if (m1 & 1) m1pm1 = -1.0;

    if (n == p) {
	return (15.0/L4*(m1pm1 + 1.0)) * (L5*(4.0*PI5*n5 + 20.0*PI3*n3 - 30.0*M_PI*dn)) / (40.0*PI5*n5);
    }
    return -(15.0/L4*(m1pm1 + 1.0)) * (L*((4.0*M_PI*L*L*cos(M_PI*npp)*(6.0*L*L - L*L*PI2*npp2))/npp4 + (4.0*M_PI*L*L*cos(M_PI*nmp)*(6.0*L*L - L*L*PI2*nmp2))/nmp4))/(2.0*PI5);
}

double RectangularImperfectPlate::i11_mat(int Npsi, int Nphi, double L, int m, int n, int p)
{
    int m1 = m - 1;
    double L3 = L*L*L;
    double L4 = L3*L;
    double dp = (double)p;
    double p2 = dp*dp;
    double nmp = (double)n - (double)p;
    double nmp2 = nmp*nmp;
    double nmp4 = nmp2*nmp2;
    double npp = (double)n + (double)p;
    double npp2 = npp*npp;
    double npp4 = npp2*npp2;
    double m1pm1 = 1.0;
    if (m1 & 1) m1pm1 = -1.0;

    if (n == p) {
	return (-4.0/L3*(7.0*m1pm1 + 8.0))*L4/8.0 + (-4.0/L3*(7.0*m1pm1 + 8.0))*(3.0*L4)/(8.0*PI2*p2);
    }
    return (-4.0/L3*(7.0*m1pm1 + 8.0))*(L*((6.0*L3)/nmp4 + (6.0*L3)/npp4))/(2.0*PI4) - (-4.0/L3*(7.0*m1pm1 + 8.0))*(L*((3.0*L*cos(M_PI*npp)*(2.0*L*L - L*L*PI2*npp2))/npp4 + (3.0*L*cos(M_PI*nmp)*(2.0*L*L - L*L*PI2*nmp2))/nmp4))/(2.0*PI4);
}

double RectangularImperfectPlate::i12_mat(int Npsi, int Nphi, double L, int m, int n, int p)
{
    int m1 = m - 1;
    double L3 = L*L*L;
    double dp = (double)p;
    double p2 = dp*dp;
    double npp = (double)n + (double)p;
    double npp2 = npp*npp;
    double nmp = (double)n - (double)p;
    double nmp2 = nmp*nmp;
    double m1pm1 = 1.0;
    if (m1 & 1) m1pm1 = -1.0;

    if (n == p) {
	return (6.0/(L*L)*(2.0*m1pm1 + 3.0))*L3/6.0 + (6.0/(L*L)*(2.0*m1pm1 + 3.0))*L3/(4.0*PI2*p2);
    }
    return (6.0/(L*L)*(2.0*m1pm1 + 3.0))*L3*cos(M_PI*nmp)/(PI2*nmp2) + (6.0/(L*L)*(2.0*m1pm1 + 3.0))*L3*cos(M_PI*npp)/(PI2*npp2);
}

double RectangularImperfectPlate::i13_mat(int Npsi, int Nphi, double L, int m, int n, int p)
{
    if (n == p) return -L/2.0;
    return 0.0;
}

int RectangularImperfectPlate::AiryStressFactorsCalculation(int Npsi, double Lx, double Ly, double *coeff0, double *coeff1, double *coeff2)
{
    int m, n, p, q, r, s;
    int Npsi2 = Npsi*Npsi;
    int S;
    int bestidx;

    double *K, *M;
    double *L, *U, *Linv, *Uinv, *C, *tmp0, *tmp1, *tmp2;
    double *d, *v;
    double *VEC, *VAL;
    double bestval = 0.0;
    double norm;

    double *NN, *MM, *nmatr, *nmatr2;

    /* allocate all temporary storage needed */
    K = new double[Npsi2*Npsi2];
    M = new double[Npsi2*Npsi2];
    L = new double[Npsi2*Npsi2];
    U = new double[Npsi2*Npsi2];
    Linv = new double[Npsi2*Npsi2];
    Uinv = new double[Npsi2*Npsi2];
    C = new double[Npsi2*Npsi2];
    tmp1 = new double[Npsi2*Npsi2];
    if (coeff0 != NULL){ 
    	tmp0 = new double[Npsi2*Npsi2];
    	tmp2 = new double[Npsi2*Npsi2];
    }
    d = new double[Npsi2];
    v = new double[Npsi2*Npsi2];
    VAL = new double[Npsi2];
    VEC = new double[Npsi2*Npsi2];

    /* calculate K and M matrices */
    r = 0;
    for (m = 0; m < Npsi; m++) {
	for (n = 0; n < Npsi; n++) {
	    s = 0;
	    for (p = 0; p < Npsi; p++) {
		for (q = 0; q < Npsi; q++) {
		    K[(r*Npsi2)+s] = int1(m,p,Lx)*int2(n,q,Ly) + int2(m,p,Lx)*int1(n,q,Ly) +
			2.0*int4(m,p,Lx)*int4(n,q,Ly);
		    M[(r*Npsi2)+s] = int2(m,p,Lx)*int2(n,q,Ly);
		    s++;
		}
	    }
	    r++;
	}
    }

    /* solve generalised Eigenvalue problem for K and M */
    /* get Cholesky decomposition of M */
    denseCholeskyDecomp(Npsi*Npsi, M, L, U);

    /* get inversion of lower triangle */
    invertLowerTriangle(Npsi*Npsi, L, Linv);

    /* get transpose of inversion */
    transposeDenseMatrix(Npsi*Npsi, Linv, Uinv);

    /* tmp1 = Linv*K */
    denseMatrixMatrixMultiply(Npsi*Npsi, Linv, K, tmp1);

    /* C = Linv*K*Uinv */
    denseMatrixMatrixMultiply(Npsi*Npsi, tmp1, Uinv, C);

    /* compute eigenvalues of C */
    getEigenvalues(Npsi*Npsi, C, d, v);


    /* sort eigenvalues into order, remove negative or zero ones. also sort eigenvectors */
    S = 0;
    do {
	bestval = 1e40;
	bestidx = -1;
	for (m = 0; m < Npsi2; m++) {
	    if ((d[m] < bestval) && (d[m] > 0.0)) {
		bestval = d[m];
		bestidx = m;
	    }
	}
	
	if (bestidx >= 0) {
	    /* found one */
	    VAL[S] = bestval;
	    d[bestidx] = 0.0;

	    /* copy the vector as well */
	    for (n = 0; n < Npsi2; n++) {
		tmp1[n] = v[(n*Npsi2) + bestidx];

	    }

	    /* multiply the vector by Uinv to correct it */
	    denseMatrixVectorMultiply(Linv, tmp1, &VEC[S*Npsi2], Npsi2, Npsi2);

	    /* normalise the vector */
	    norm = 0.0;
	    for (n = 0; n < Npsi2; n++) {
		norm += (VEC[S*Npsi2+n]*VEC[S*Npsi2+n]);	
	    }
	    norm = sqrt(norm);
	    for (n = 0; n < Npsi2; n++) {
		VEC[S*Npsi2+n] /= norm;
	    }

	    S++;
	}
    } while (bestidx >= 0);

    /* create the two integral matrices */
    NN = int2_mat(Npsi, Lx); /* these are both symmetric */
    MM = int2_mat(Npsi, Ly);

    /* NN is a single column, then replicated Npsi2 times */
    /* MM is a single row, then replicated Npsi2 times */
    /* element-wise multiply them to get nmatr in *column-major order* */
    nmatr = new double[Npsi2*Npsi2];
    for (m = 0; m < Npsi2; m++) {
	for (n = 0; n < Npsi2; n++) {
	    nmatr[(n*Npsi2)+m] = NN[m] * MM[n];
	}
    }
    
    /* now we treat nmatr as a 4D NpsixNpsixNpsixNpsi array and permute the Nphiensions into the order:
     * 4 1 3 2 */
    nmatr2 = new double[Npsi2*Npsi2];
    for (m = 0; m < Npsi; m++) {
	for (n = 0; n < Npsi; n++) {
	    for (p = 0; p < Npsi; p++) {
		for (q = 0; q < Npsi; q++) {
		    nmatr2[q + (m*Npsi) + (p*Npsi2) + (n*Npsi2*Npsi)] =
			nmatr[m + (n*Npsi) + (p*Npsi2) + (q*Npsi2*Npsi)];
		}
	    }
	}
    }

    /* from here on nmatr(2) is treated as a row vector */

    /* loop over eigenvectors */
    for (m = 0; m < S; m++) {
	double svm = sqrt(VAL[m]);
	double temp, temp2, temp3, norm, snorm;

	norm = 0.0;
	for (p = 0; p < Npsi2; p++) { /* columns of temp3 */
	    for (q = 0; q < Npsi2; q++) { /* rows of temp3 */
		temp = VEC[(m*Npsi2)+q];
		temp2 = VEC[(m*Npsi2)+p];
		temp3 = temp * temp2;
		norm += temp3 * nmatr2[(p*Npsi2)+q];
	    }
	}

	snorm = sqrt(norm);
	for (n = 0; n < Npsi2; n++) {		    	
	    tmp1[(m*Npsi2)+n] = VEC[(m*Npsi2)+n] / snorm / svm;	    
	    if (coeff0 != NULL){
	    	tmp0[(m*Npsi2)+n] = VEC[(m*Npsi2)+n] / snorm;
	    	tmp2[(m*Npsi2)+n] = VEC[(m*Npsi2)+n] / snorm / svm / svm;
	    }
	}
    }

    /* compute actual S value */
    S = S / 2;

    /* copy result into actual buffer */
    for (m = 0; m < S; m++) {
	for (n = 0; n < S; n++) {
		coeff1[(n*S)+m] = tmp1[(n*Npsi2)+m];
		if (coeff0 != NULL){
			coeff0[(n*S)+m] = tmp0[(n*Npsi2)+m];
	    	coeff2[(n*S)+m] = tmp2[(n*Npsi2)+m];
		}
	}
    }

    delete[] K;
    delete[] M;
    delete[] L;
    delete[] U;
    delete[] Linv;
    delete[] Uinv;
    delete[] C;
    delete[] tmp0;
    delete[] tmp1;
    delete[] tmp2;
    delete[] d;
    delete[] v;
    delete[] VAL;
    delete[] VEC;
    delete[] NN;
    delete[] MM;
    delete[] nmatr;
    delete[] nmatr2;

    return S;
}


/* integration functions */
double RectangularImperfectPlate::int1(int m, int p, double L)
{
    double y;
    if ((m == 0) && (p == 0)) {
	y = 720.0 / (L*L*L);
    }
    else if (m == p) {
	double m1pm = 1.0;
	if (m & 1) m1pm = -1.0; /* (-1)^m */

	y = (PI4 * ((double)m*m*m*m) - 672.0*m1pm - 768.0) / (2.0*L*L*L);
    }
    else if ((m == 0) || (p == 0)) {
	y = 0.0;
    }
    else {
	double m1pm = 1.0;
	double m1pp = 1.0;
	if (m & 1) m1pm = -1.0; /* (-1)^m */
	if (p & 1) m1pp = -1.0; /* (-1)^p */

	y = -(24.0*(7.0*m1pm + 7.0*m1pp + 8.0*m1pm*m1pp + 8.0)) / (L*L*L);
    }
    return y;
}

double RectangularImperfectPlate::int2(int m, int p, double L)
{
    double y;
    double m4, p4;
    double m1pm = 1.0;
    double m1pp = 1.0;
    if (m & 1) m1pm = -1.0; /* (-1)^m */
    if (p & 1) m1pp = -1.0; /* (-1)^p */

    m4 = (double)m;
    m4 = m4*m4*m4*m4;
    p4 = (double)p;
    p4 = p4*p4*p4*p4;

    if ((m == 0) && (p == 0)) {
	y = (10.0*L)/7.0;
    }
    else if (m == p) {
	y = (67.0*L)/70.0 - (m1pm*L)/35.0 - (768.0*L)/(PI4*m4) - (672.0*m1pm*L)/(PI4*m4);
    }
    else if (m == 0) {
	y = (3.0*L*(m1pp+1.0)*(PI4*p4 - 1680.0)) / (14.0*PI4*p4);
    }
    else if (p == 0) {
	y = (3.0*L*(m1pm+1.0)*(PI4*m4 - 1680.0)) / (14.0*PI4*m4);
    }
    else {
	y = -(L*(11760.0*m1pm + 11760.0*m1pp - 16.0*PI4*m4 + 13440.0*m1pm*m1pp +
		 m1pm*PI4*m4 + m1pp*PI4*m4 - 16.0*m1pm*m1pp*PI4*m4 + 13440.0)) / (70.0*PI4*m4) -
	    (L*(13440.0*m4 + 11760.0*m1pm*m4 + 11760.0*m1pp*m4 + 13440.0*m1pm*m1pp*m4)) / (70.0*PI4*m4*p4);
    }
    return y;
}

double RectangularImperfectPlate::int4(int m, int p, double L)
{
    double y;
    double m2 = (double)(m*m);
    double p2 = (double)(p*p);
    double m1pm = 1.0;
    double m1pp = 1.0;
    if (m & 1) m1pm = -1.0; /* (-1)^m */
    if (p & 1) m1pp = -1.0; /* (-1)^p */

    if ((m == 0) && (p == 0)) {
	y = 120.0 / (7.0*L);
    }
    else if ((m == p) && (m != 0)) {
	y = (768.0*PI2*m2 - 47040.0*m1pm + 35.0*PI4*m2*m2 + 432.0*m1pm*PI2*m2 - 53760.0) /
	    (70.0*L*PI2*m2);
    }
    else if (m == 0) {
	y = (60.0*(m1pp + 1.0)*(PI2*p2 - 42.0)) / (7.0*L*PI2*p2);
    }
    else if (p == 0) {
	y = (60.0*(m1pm + 1.0)*(PI2*m2 - 42.0)) / (7.0*L*PI2*m2);
    }
    else {
	y = 192.0/35.0/L*(1.0 + m1pm*m1pp) - 192.0/m2/p2/L/PI2*((p2+m2)*(1.0+m1pm*m1pp)) -
	    168.0/m2/p2/L/PI2*((p2+m2)*(m1pm+m1pp)) + 108.0/35.0/L*(m1pm+m1pp);
    }
    return y;
}

double *RectangularImperfectPlate::int2_mat(int tt, double L)
{
    double *y;
    int m, p;
    double m1pm, m1pp, m4;
    double dp;

    y = new double[tt*tt];

    for (m = 0; m < tt; m++) {
	m4 = (double)m;
	m4 = m4*m4*m4*m4;
	m1pm = 1.0;
	if (m & 1) m1pm = -1.0;

	for (p = 1; p < tt; p++) {
	    dp = (double)p;
	    m1pp = 1.0;
	    if (p & 1) m1pp = -1.0;

	    if (m == 0) {
		/* top row */
		y[(m*tt)+p] = (3.0*L*(m1pp + 1.0) * (PI4*(dp*dp*dp*dp) - 1680.0)) / (14.0*PI4*(dp*dp*dp*dp));
	    }
	    else {
		y[(m*tt)+p] = -(L*(11760.0*m1pm + 11760.0*m1pp - 16.0*PI4*m4 +
				   13440.0*m1pm*m1pp + m1pm*PI4*m4 + m1pp*PI4*m4 -
				   16.0*m1pp*PI4*(m4*m1pm) + 13440.0))/(70.0*PI4*m4) -
		    (L*(13440.0*m4 + 11760.0*m1pm*m4 + 11760*m1pp*m4 + 13440.0*m1pm*m1pp*m4)) /
		    (70.0*PI4*m4*dp*dp*dp*dp);
	    }
	}

	/* diagonal element */
	y[(m*tt)+m] = (67.0*L)/70.0 - (m1pm*L)/35.0 - (768.0*L)/(PI4*m4) -
	    (672.0*m1pm*L)/(PI4*m4);

	/* left hand side */
	y[(m*tt)] = (3.0*L*(m1pm + 1.0) * (PI4*m4 - 1680.0)) / (14.0*PI4*m4);
    }

    /* top left */
    y[0] = (10.0*L)/7.0;

    return y;
}


void RectangularImperfectPlate::LoadHTensor(double *coeff0, double *coeff1, double *coeff2, int Nphi, int Npsi, double Lx,
		double Ly, int S, double *H0, double *H1, double *H2) 
{
	char filenameH0[1000], filenameH1[1000], filenameH2[1000];		
	int NphiF, NpsiF;
	bool HLoaded = false;

	//double e = (Lx*Ly / 4.0) / rho * Young ;
	double Hscale = sqrt(e / 2.0);
		
	GlobalSettings *gs = GlobalSettings::getInstance();
	
	sprintf(filenameH0, "./RectangularPlateData/H0_SS-Nphi_%d-S_%d-Lx_%f-Ly_%f.bin", Nphi, S, Lx, Ly);
	sprintf(filenameH1, "./RectangularPlateData/H1_SS-Nphi_%d-S_%d-Lx_%f-Ly_%f.bin", Nphi, S, Lx, Ly);
	sprintf(filenameH2, "./RectangularPlateData/H2_SS-Nphi_%d-S_%d-Lx_%f-Ly_%f.bin", Nphi, S, Lx, Ly);
	

	
	
	while (!HLoaded){
		
		if (gs->getStormerVerletScheme()){
			ifstream fiH1(filenameH1, ios::in | ios::binary);
			if (fiH1.good()) {
					logMessage(1, "Loading H1 from file");
				
					
					double *H1_aux = new double;
					int p;
					
						for (int i=0; i<Npsi; i++){
							p = 0;
							for (int j=0; j<Nphi*Nphi; j++){										
								fiH1.seekg((i+j*S)*sizeof(double), ios::beg);
								fiH1.read((char *)H1_aux,sizeof(double));
								H1[i+p*Npsi] = *H1_aux * Hscale;
								p++; 	
							}
						}
						HLoaded = true;
						fiH1.close();

						delete H1_aux;	


			} 
			else {

				
				H_tensorRectangular(coeff0, coeff1, coeff2, Nphi, Npsi, Lx, Ly, S, H0, H1, H2);	
				
				
				
				logMessage(1, "Performed H_tensorRectangular");
				
				
			}
			} 
			else {
				ifstream fiH0(filenameH0, ios::in | ios::binary);
				ifstream fiH1(filenameH1, ios::in | ios::binary);
				ifstream fiH2(filenameH2, ios::in | ios::binary);
				
				if (fiH0.good() && fiH1.good() && fiH2.good()) {
						logMessage(1, "Loading H tensors from file");
						
						
						double *H0_aux = new double;
						double *H1_aux = new double;
						double *H2_aux = new double;
						int p;
						
							for (int i=0; i<Npsi; i++){
								p = 0;
								for (int j=0; j<Nphi*Nphi; j++){										
									fiH0.seekg((i+j*S)*sizeof(double), ios::beg);
									fiH0.read((char *)H0_aux,sizeof(double));
									H0[i+p*Npsi] = *H0_aux * Hscale;
									
									fiH1.seekg((i+j*S)*sizeof(double), ios::beg);
									fiH1.read((char *)H1_aux,sizeof(double));
									H1[i+p*Npsi] = *H1_aux * Hscale;
									
									fiH2.seekg((i+j*S)*sizeof(double), ios::beg);
									fiH2.read((char *)H2_aux,sizeof(double));
									H2[i+p*Npsi] = *H2_aux * Hscale;
									
									p++; 	
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
					H_tensorRectangular(coeff0, coeff1, coeff2, Nphi, Npsi, Lx, Ly, S, H0, H1, H2);	
					
					logMessage(1, "Performed H_tensorRectangular");
					
					
				}
			}
	}
		
		
		
}

void RectangularImperfectPlate::RectangularImperfection(double *Imperfection, double Lx, double Ly, double H, int Nx, int Ny, char ImperfectionType, double xWidth, double yWidth)
{
	int x, y;

	
	double Y0, X0;
	
	double xVec, yVec, fsx, fsy ;
	
	
	int nx_min, nx_max, ny_min, ny_max;

		
	memset(Imperfection, 0, Nx*Ny* sizeof(double));

	switch (ImperfectionType){
		case 'r': // 2D Raised Cosine
			// x direction

			X0 = round(Nx/2);
			fsx = (double)Nx/Lx;
			nx_min = floor(X0 - fsx*xWidth);
			nx_max = floor(X0 + fsx*xWidth);

			// y direction
			Y0 = round(Ny/2);
			fsy = (double)Ny/Ly;
			ny_min = floor(Y0 - fsy*yWidth);
			ny_max = floor(Y0 + fsy*yWidth);


			for (int nx = nx_min; nx<nx_max; nx++){
				xVec = 0.5*(1 + cos(M_PI*(((double)nx-1)/fsx/xWidth - X0/fsx/xWidth)));
				
				for (int ny = ny_min; ny<ny_max; ny++){
					yVec = 0.5*(1 + cos(M_PI*(((double)ny-1)/fsy/yWidth - Y0/fsy/yWidth)));
					
					Imperfection[nx*Ny + ny] = H*xVec*yVec;

				}
				
			}
			
		break;
	}
	
	
	
}


void RectangularImperfectPlate::ProjectionCoefficients(double *Imperfection, double *Ai, int Nphi, int modeType, double max_error, int Nx, int Ny, double Lx, double Ly)
{

	double zg;
	bool Include;
	double err_i, err_std, err_av;
	int i=0;
	int j;
	double norm=0;

	double *recons = new double[Nx*Ny];
	double *phi = new double[Nx*Ny];

	
	zg = 0; // 0 for axisymmetric caps

	memset(recons, zg, Nx*Ny * sizeof(double));
		
	err_i = 99999;
	

	
	while ((err_i>max_error)&&(i<Nphi)){

		Include = false;


		switch (modeType){
			case 0:
				Include = true;
				break;			
		}

		if (Include == true){

			ModeShape(phi, BC, kx[i], ky[i], Lx, Ly, Nx, Ny);
		
			norm = ScalarProductCartesian(phi,phi, Lx, Ly, Nx, Ny);

			
			Ai[i] = ScalarProductCartesian(Imperfection, phi, Lx, Ly, Nx, Ny)/sqrt(norm);

			
			err_std = 0;
			err_av = 0;
			for (j=0; j<(Nx*Ny); j++){
				recons[j] = recons[j] + Ai[i]*phi[j]/sqrt(norm); //phi should be normalized before computing Ai. Instead, it is divided by norm twice here.
				err_av = err_av + abs((Imperfection[j]-recons[j])/Imperfection[j]);		

			}
			
			err_av = err_av / (Nx*Ny);
			
			if ((Ai[i] != 0) || (err_av != 1 )){
				for (j=0; j<(Nx*Ny); j++){
					err_std = err_std + (abs((Imperfection[j]-recons[j])/Imperfection[j]) - err_av)*(abs((Imperfection[j]-recons[j])/Imperfection[j]) - err_av);	
				}
				
				err_std = sqrt(err_std / (Nx*Ny - 1));
				
				if (err_i > err_std){
					err_i = err_std;
				}
			}

		} else{
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

void RectangularImperfectPlate::ModeShape(double *phi, char BC, int kx, int ky, double Lx, double Ly, int Nx, int Ny)
{
	
	double dx, dy, x, y;

	
	dx = Lx / (Nx - 1);
	dy = Ly / (Ny - 1);
	
	if (BC == 's'){
	for (int nx = 0; nx<Nx; nx++ ){
		x = dx*nx;
		for (int ny = 0; ny<Ny; ny++ ){
			y = dy*ny;
			phi[nx*Ny + ny] = sin(kx*M_PI*x/Lx)*sin(ky*M_PI*y/Ly);

			
		}
	}
	}


}

