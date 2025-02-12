/*
 * VK-Gong code 
 *
 * Author of this version: Àngels Aragonès
 *
 * Based on VK-Gong for Matlab and NeSS Framework Code (Copyright (c) The University of Edinburgh, 2014-2015. All rights reserved. Author: James Perry)
 *
 * Non-linear modal imperfect circular plate component
 */

#ifndef _CIRCULAR_IMPERFECT_PLATE_H_
#define _CIRCULAR_IMPERFECT_PLATE_H_

#include "Component.h"




// For the moment, this class inherits from Component. This may change to ModalPlate.
class CircularImperfectPlate : public Component {
 public:
	CircularImperfectPlate(string name, int iNphi, int iNpsi, double Rd, double h, double H, char ImperfectionType, double tau2, int modeType,
			       double nu, double Young, double rho, char BC, double KR, double KT, int Nr, int Nth, double dFac, double dExp, double dCons, int fs);
    virtual ~CircularImperfectPlate();
    
    double *getOmega() { return ov; }
    double *getXi() { return xi;}
    int *getk() {return kt;}
    int *getn() {return nt;}
    int *getc() {return ct;}
    double geth() { return h; }
    double getRho() { return rho; }
    double getHeight() { return H; }
    double getRd() { return Rd; }
    double getnu() { return nu; }
    int getNr() {return Nr;}
    int getNth() {return Nth;}
    double *getCnorm() { return Cnorm; }
    double gettnd() { return tnd; }
    double getYoung() { return Young;}
    double gete() { return e;}
    double getKR() { return KR;}
    char getBoundaryConditions() { return BC;}
    double Norm_Modes(int k, double xkn, double R);


	typedef struct ModeXKN{
			double x;
			int k;
			int n;
		    bool operator < (const ModeXKN& rhs) const {
		        return x < rhs.x;
		    }
	}ModeXKN;
	
	
    
    
 protected:
    void plate_def(int *fs);
    void TransverseModeCalculator(double dx, double xmax, char BC, double nu, double KR, double KT, int Nphi, double* xi, int* kt, int* ct, double* ov);
    void InPlaneModeCalculator(double dx, double xmax, char BC, double nu, int Npsi, double* zl, int* kl, int* cl);
    
	double* SortZeros( int xkLength, int Total, std::vector<CircularImperfectPlate::ModeXKN> xk, char BC);

	void LoadHTensor(double dr_H);
	void H_tensorCircular(int NphiF, int NpsiF, char BC, double nu, double KR, double KT, double dr_H );
	double HCoefficientCircular( int kp, int kq, int cp, int cq, double xip, double xiq, int ki, int ci, double zeta, double nu, double KR,  double dr_H );
	
    void AxisymmetricCap(double H, double Rd, char ImperfectionType, double *Imperfection, int Nr, int Nth);
    void ProjectionCoefficients(double *Imperfection, double *Ai, int Nphi, int modeType, double max_error,  int Nr, int Nth, double nu, double KR, int *k, int *n, int *c, double *xkn);
    void ModeShape(double *phi, int k,int c, double xkn, double Rd, double nu, double KR, int Nr, int Nth);



    double CosCosCosIntegration(int k, int l, int m);
    double CosSinSinIntegration(int k, int l, int m);

    // dimensions in metres
    double h;
    double H;
    double Rd;
    double nu;
    double Young;
    double rho;
    
    double dFac;
    double dExp;
	double dCons;

    char ImperfectionType;
    double tau2;
    int modeType;
    int Nr;
    int Nth;


    double *ov, *xi; // Dimensionless omega and zeta of the transverse modes where zeta^4=omega^2
    int *kt, *nt, *ct; // Nodal radius, nodal circles and cos/sin shape of transverse modes.

    double *zl; // Dimensionless omega and zeta of the in-plane modes where zeta^4=omega^2
    int *kl, *nl, *cl; // Nodal radius, nodal circles and cos/sin shape of in-plane modes.

    double *Cnorm;

    int *ModeIndices;
    
    double e;
    double D;
    double tnd;
    
    char BC;
    double KR;
    double KT;


};

#endif
