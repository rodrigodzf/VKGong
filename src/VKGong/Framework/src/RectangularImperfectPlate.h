/*
 * VK-Gong code 
 *
 * Author of this version: Àngels Aragonès
 *
 * Based on VK-Gong for Matlab and NeSS Framework Code (Copyright (c) The University of Edinburgh, 2014-2015. All rights reserved. Author: James Perry)
 *
 * Non-linear modal imperfect rectangular plate component
 */

#ifndef _RECTANGULAR_IMPERFECT_PLATE_H_
#define _RECTANGULAR_IMPERFECT_PLATE_H_

#include "Component.h"


class RectangularImperfectPlate : public Component {
 public:
    RectangularImperfectPlate(string name, int iNphi, int iNpsi, double Lx, double Ly, double h, double H, char ImperfectionType,
		       double xWidth, double yWidth, int modeType, double nu, double Young, double rho, char BC, int Nx, int Ny, double dFac, double dExp, double dCons, int fs);
    virtual ~RectangularImperfectPlate();

   
    double getLx() { return Lx; }
    double getLy() { return Ly; }
    double *getOmega() { return ov; }
    int *getKx() { return kx; }
    int *getKy() { return ky; }
    double geth() { return h; }
    double getRho() { return rho; }
    double *getCnorm() { return Cnorm; }
    double gete() { return e; }
    char getBoundaryConditions() { return BC; }


 protected:
    void plate_def(int *fs);

    int AiryStressFactorsCalculation(int A, double Lx, double Ly, double *coeff0, double *coeff1, double *coeff2);
    double int1(int m, int p, double L);
    double int2(int m, int p, double L);
    double int4(int m, int p, double L);
    double *int2_mat(int tt, double L);

    void H_tensorRectangular(double *coeff0, double *coeff1, double *coeff2, int DIM, int A, double Lx, double Ly,
			int S, double *H0, double *H1, double *H2);
    void LoadHTensor(double *coeff0, double *coeff1, double *coeff2, int Nphi, int Npsi, double Lx,
    		double Ly, int S, double *H0, double *H1, double *H2);
    void ComputeTransverseEigenfrequenciesRectangular(int Nphi, double Lx, double Ly, double h, double Young, double rho, double nu, double* ov, int* kx, int* ky);
    
    void RectangularImperfection(double *Imperfection, double Lx, double Ly, double H, int Nx, int Ny, char ImperfectionType, double xWidth, double yWidth);
    void ProjectionCoefficients(double *Imperfection, double *Ai, int Nphi, int modeType, double max_error, int Nx, int Ny);
    void ProjectionCoefficients(double *Imperfection, double *Ai, int Nphi, int modeType, double max_error, int Nx, int Ny, double Lx, double Ly);
    void ModeShape(double *phi, char BC, int kx, int ky, double Lx, double Ly, int Nx, int Ny);
    
    double i1_mat(int A, int DIM, double L, int m, int n, int p);
    double i2_mat(int A, int DIM, double L, int m, int n, int p);
    double i3_mat(int A, int DIM, double L, int m, int n, int p);
    double i4_mat(int A, int DIM, double L, int m, int n, int p);
    double i5_mat(int A, int DIM, double L, int m, int n, int p);
    double i9_mat(int A, int DIM, double L, int m, int n, int p);
    double i10_mat(int A, int DIM, double L, int m, int n, int p);
    double i11_mat(int A, int DIM, double L, int m, int n, int p);
    double i12_mat(int A, int DIM, double L, int m, int n, int p);
    double i13_mat(int A, int DIM, double L, int m, int n, int p);

    // dimensions in metres
    double Lx, Ly;


    double h;
    double nu;
    double Young;
    double rho;

    double *ov;
    int *kx, *ky;
    double *Cnorm;
    
    double H;
    
    double dFac, dExp, dCons;
    char ImperfectionType, BC;
    int modeType;
    
    double e;
    
    int Nx;
    int Ny;
    
    float xWidth, yWidth;


};

#endif
