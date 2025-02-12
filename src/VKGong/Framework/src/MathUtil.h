/*
 * VK-Gong code 
 *
 * Author of this version: Àngels Aragonès
 *
 * Based on VK-Gong for Matlab and NeSS Framework Code (Copyright (c) The University of Edinburgh, 2014-2015. All rights reserved. Author: James Perry)
 *
 * Useful maths functions
 */

#ifndef _MATH_UTIL_H_
#define _MATH_UTIL_H_



extern "C" {
#include "csrmatrix.h"
};

/*
 * Holds all the settings and temporary buffers used by the vector Newton solver
 */
typedef struct newton_solver_s
{
    int maxIterations;
    int jacobiIterations;
    double tolerance;
    double jacobiTolerance;

    CSRmatrix *IMQ;
    CSRmatrix *J;

    double *coeff;
    double *phi_a;
    double *fac2;
    double *fac3;
    double *temp;
    double *F;
    double *IMQdiag;
    double *x;

    double *d;
    int *jacobiRows;
    double *xp;
} newton_solver_t;

/*
 * Create and initialise a vector Newton solver structure. The matrix M must have the same
 * non-zero structure as the matrix that is later passed to newtonSolverVector (with no diagonal),
 * but it doesn't have to contain the same values.
 */
newton_solver_t *newtonSolverCreate(int max_it, double tol, int jacobi_it, double jacobi_tol,
				    CSRmatrix *M);

/* Free vector Newton solver structure */
void newtonSolverFree(newton_solver_t *newton);

/* vector Newton solver */
void newtonSolverVector(newton_solver_t *newton, double *r, double *a, double *b, CSRmatrix *M, double *Mdiag,
			double *q, double *K, double *alpha, double *phi_ra, double *R);

/* scalar Newton solver */
double newtonSolver(double a, double b, double M, double K, double alpha, double offset, double one_sided,
		    double *phi_ra);

#endif

/* implementation of Matlab's interp1 function */
void interp1(double *x, double *v, double *xq, int lenx, int lenxq, double *result);

/* as above, but x and v are interleaved in a single array */
void interp1_interleaved(double *xandv, double *xq, int lenx, int lenxq, double *result);

double factorial(double n);

/*
 * Computes the Crout LU decomposition of a square dense matrix A (size n)
 */
bool croutDecomposition(double *A, double *L, double *U, int n);

/*
 * Performs a triangular solve using the Crout decomposition
 * y is temporary storage
 */
void croutSolve(double *L, double *U, double *b, double *x, double *y, int n);

/*
 * Gets the eigenvalues (in val) and the corresponding eigenvectors (in vec)
 * of a dense symmetric matrix A, size NxN, using the Jacobi method.
 * The eigenvalues are not sorted by this function.
 * A is overwritten in the process.
 */
void getEigenvalues(int N, double *A, double *val, double *vec);

/*
 * Multiplies m-by-n matrix A by vector x, stores result in b
 */
void denseMatrixVectorMultiply(double *A, double *x, double *b, int m, int n);

/*
 * As above, but the matrix is transposed first
 * The dimensions given are the non-transposed dimensions
 */
void denseMatrixVectorMultiplyTransposed(double *A, double *x, double *b, int m, int n);

/*
 * Performs Cholesky decomposition of B into L and U. All matrices are NxN
 */
void denseCholeskyDecomp(int N, double *B, double *L, double *U);

/* inverts lower triangular matrix L into I */
void invertLowerTriangle(int N, double *L, double *I);

/* transposes dense NxN matrix in into out */
void transposeDenseMatrix(int N, double *in, double *out);

/* transposes dense NrxNc matrix in into out */
void transposeMatrixRectangular(int Nr, int Nc, double *in, double *out);

/* out = A * B; */
void denseMatrixMatrixMultiply(int N, double *A, double *B, double *out);
void denseMatrixMatrixMultiplyRectangular(int Nr, int Nc, double *A, double *B, double *out);


double Trapz2Dpolar(double* f , double r_min, double r_max, double th_min, double th_max, int Nr, int Nth );
double Trapz1DPolar(double* f , double r_min, double r_max, int Nr );
double ScalarProductPolar(double *f1, double *f2, int Nr, int Nth);

double Trapz2DCartesian(double* f , double x_min, double x_max, double y_min, double y_max, int Nx, int Ny );
double Trapz1DCartesian(double* f , double x_min, double x_max, int pointsx );
double ScalarProductCartesian(double *f1, double *f2, double Lx, double Ly, int Nx, int Ny);

//bool gaussSolver(double* mat, double* rhs, double* minPivot, int* row, int nbRows);
bool gaussSolver(double* mat, double* rhs, int nbRows);
