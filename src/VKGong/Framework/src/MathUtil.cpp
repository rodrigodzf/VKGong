/*
 * VK-Gong code 
 *
 * Author of this version: Àngels Aragonès
 *
 * Based on VK-Gong for Matlab and NeSS Framework Code (Copyright (c) The University of Edinburgh, 2014-2015. All rights reserved. Author: James Perry)
 *
 */

#include "MathUtil.h"

#include <cmath>
#include <cstring>
using namespace std;


#define EPSILON 2.2204460492503131e-16

static double sign(double val)
{
    if (val < 0.0) return -1.0;
    else if (val > 0.0) return 1.0;
    else return 0.0;
}

double newtonSolver(double a, double b, double M, double K, double alpha, double offset, double one_sided,
		    double *phi_ra)
{
    double coeff = K / (alpha + 1.0);
    double r = 1.0;
    double R, F, temp;
    double phi_ra2;

    for (int nn = 0; nn < 5; nn++) {
	double ra = r + a;
	double rae = fabs(ra) - offset;
	double ae = fabs(a) - offset;

	double sra = sign(ra);
	double srae = sign(rae);
	double sae = sign(ae);
	double sa = sign(a);
	double sr = sign(r);

	phi_ra2 = 0.5 * coeff * (1.0 - 0.5 * one_sided * (1.0 - sra)) * (1.0 + srae)
	    * pow(fabs(rae), alpha + 1.0);
	double phi_a = 0.5 * coeff * (1.0 - 0.5 * one_sided * (1.0 - sa)) * (1.0 + sae)
	    * pow(fabs(ae), alpha + 1.0);
	double phi_prime = 0.5 * K * sra * (1.0 - 0.5 * one_sided * (1.0 - sra))
	    * (1.0 + srae) * pow(fabs(rae), alpha);
	if (fabs(r) > EPSILON) {
	    R = sr * (phi_ra2 - phi_a) / fabs(r);
	    F = r + M*R + b;
	    temp = (r * phi_prime - phi_ra2 + phi_a) / (r*r);
	}
	else {
	    R = 0.5 * (K * (1.0 - 0.5 * one_sided * (1.0 - sae))) * (1.0 + sae)
		* pow(fabs(ae), alpha);
	    F = r + M*R + b;
	    temp = 0.5 * (0.5 * alpha * (alpha+1.0) * K
			  * (1.0 - 0.5 * one_sided * (1.0 - sae)) * (1.0 + sae)
			  * pow(fabs(ae), (alpha - 1.0)));
	}
	r = r - F / (1.0 + M * temp);
	//printf("%f %f %f %f %f %f %f %f %f %f\n", ra, rae, ae, phi_ra2, phi_a, phi_prime,
	//       R, F, temp, r);
    }
    if (phi_ra) *phi_ra = phi_ra2;
    return -R;
}

/*
 * Create and initialise a Newton solver structure. The matrix M must have the same
 * non-zero structure as the matrix that is later passed to newtonSolver (with no diagonal),
 * but it doesn't have to contain the same values.
 */
newton_solver_t *newtonSolverCreate(int max_it, double tol, int jacobi_it, double jacobi_tol,
				    CSRmatrix *M)
{
    int len = M->nrow;
    newton_solver_t *newton = new newton_solver_t;
    if (!newton) return NULL;

    newton->maxIterations = max_it;
    newton->jacobiIterations = jacobi_it;
    newton->tolerance = tol;
    newton->jacobiTolerance = jacobi_tol;

    newton->IMQ = CSR_duplicate(M);
    newton->J = CSR_duplicate(M);

    newton->coeff = new double[len];
    newton->phi_a = new double[len];
    newton->fac2 = new double[len];
    newton->fac3 = new double[len];
    newton->temp = new double[len];
    newton->F = new double[len];
    newton->IMQdiag = new double[len];
    newton->x = new double[len];
    newton->d = new double[len];
    newton->xp = new double[len];
    newton->jacobiRows = new int[len];

    return newton;
}

void newtonSolverFree(newton_solver_t *newton)
{
    CSR_free(newton->IMQ);
    CSR_free(newton->J);
    delete(newton->coeff);
    delete(newton->phi_a);
    delete(newton->fac2);
    delete(newton->fac3);
    delete(newton->temp);
    delete(newton->F);
    delete(newton->IMQdiag);
    delete(newton->x);
    delete(newton->d);
    delete(newton->xp);
    delete(newton->jacobiRows);
    delete(newton);
}


/*
 * All parameters are inputs except:
 *  - phi_ra and R are outputs
 *  - r is both an input and an output
 *  - newton structure contains all the temporary vectors required by Newton
 *
 * The M matrix should not contain any diagonal entries; the diagonal should be in the Mdiag vector instead.
 */
void newtonSolverVector(newton_solver_t *newton, double *r, double *a, double *b, CSRmatrix *M, double *Mdiag,
			double *q, double *K, double *alpha, double *phi_ra, double *R)
{
    int len = M->nrow;
    int i, j, nn;
    double val;
    int ci;

    int k;
    double resid, rtol = newton->jacobiTolerance * newton->jacobiTolerance;
    int row, hasnz;
    double Jdiag;
    double Fval;
    int numjacobirows;

    /* compute phi_a, fac2 and fac3 */
    for (i = 0; i < len; i++) {
	newton->coeff[i] = K[i] / (alpha[i] + 1.0);

	if (a[i] > 0) {
	    double pow_a_alpham1 = pow(a[i], alpha[i]-1.0);
	    double pow_a_alpha = pow_a_alpham1 * a[i];
	    double pow_a_alphap1 = pow_a_alpha * a[i];

	    newton->phi_a[i] = newton->coeff[i] * pow_a_alphap1;
	    newton->fac2[i] = 0.5 * alpha[i] * (alpha[i]+1.0) * K[i] * pow_a_alpham1;
	    newton->fac3[i] = K[i] * pow_a_alpha;
	}
	else {
	    newton->phi_a[i] = 0.0;
	    newton->fac2[i] = 0.0;
	    newton->fac3[i] = 0.0;
	}
    }

    /*
     * compute IMQ (= I + M*Q)
     * Like M, the diagonal is stored separately
     */
    for (i = 0; i < len; i++) {
	// compute diagonal entry, including adding identity matrix
	newton->IMQdiag[i] = (Mdiag[i] * q[i]) + 1.0;

	// compute off-diagonal entries
	for (j = newton->IMQ->rowStart[i]; j < newton->IMQ->rowStart[i+1]; j++) {
	    newton->IMQ->values[j] = M->values[j] * q[newton->IMQ->colIndex[j]];
	}
    }

    /* Main Newton iteration loop */
    for (nn = 0; nn < newton->maxIterations; nn++) {

	/* compute phi_ra, R and temp */
	for (i = 0; i < len; i++) {
	    double phi_prime = 0.0, phi_diff;
	    double ra = r[i] + a[i];
	    phi_ra[i] = 0.0;
	    if (ra > 0.0) {
		double pow_ra_alpha = pow(ra, alpha[i]);
		double pow_ra_alphap1 = pow_ra_alpha * ra;
		phi_ra[i] = newton->coeff[i] * pow_ra_alphap1;
		phi_prime = K[i] * pow_ra_alpha;
	    }
	    phi_diff = phi_ra[i] - newton->phi_a[i];

	    if ((r[i] > EPSILON) || (r[i] < -EPSILON)) {
		R[i] = (phi_diff / r[i]);
		newton->temp[i] = ((r[i] * phi_prime) - phi_diff) / (r[i]*r[i]);
	    }
	    else {
		R[i] = newton->fac3[i];
		newton->temp[i] = newton->fac2[i];
	    }	    
	}

	/*
	 * This loop sets up for the Jacobi solver, performing several different operations in one go:
	 *   - computes RHS for Jacobi solver (F = IMQ*r + M*R + b)
	 *   - computes system matrix for Jacobi solver (J = I + M*(temp+q))
	 *   - computes initial residual, if small enough don't need to enter Jacobi loop at all
	 *   - computes inverse of diagonal for Jacobi solver
	 *   - makes a list of rows with non-zero off-diagonal entries for the Jacobi solver
	 * we take advantage of the fact that J, M and IMQ have identical structure...
	 */
        numjacobirows = 0;
	resid = 0.0;
	for (i = 0; i < len; i++) {
	    /* compute diagonal part of F[i] */
	    val = (newton->IMQdiag[i] * r[i]) + (Mdiag[i] * R[i]);

	    /* compute J's diagonal */
	    Jdiag = (Mdiag[i] * (newton->temp[i] + q[i])) + 1.0;

	    /* invert J's diagonal for Jacobi */
	    newton->d[i] = 1.0 / Jdiag;
	    hasnz = 0;

	    /* now handle off-diagonals */
	    for (j = M->rowStart[i]; j < M->rowStart[i+1]; j++) {
		/* non-diagonal part of F[i] */
		ci = M->colIndex[j];
		val += (M->values[j] * R[ci]) + (newton->IMQ->values[j] * r[ci]);

		/* compute J's non-diagonal values */
		newton->J->values[j] = M->values[j] * (newton->temp[ci] + q[ci]);

		/* mark this row for Jacobi */
		if (newton->J->values[j] != 0.0) hasnz = 1;
	    }

	    /* compute final value of F[i] */
	    Fval = val + b[i];

	    /* update dot product of RHS (initial residual) */
	    resid += Fval * Fval;

	    /* initialise this row for Jacobi */
	    if (hasnz) {
		/* row has non-zero off-diagonal entries, add it to the list */
		newton->jacobiRows[numjacobirows] = i;
		numjacobirows++;
		newton->x[i] = 0.0;
	    }
	    else {
		/* row only has diagonal. solve it right now */
		newton->x[i] = Fval * newton->d[i];
		newton->xp[i] = Fval * newton->d[i];
	    }

	    /* store F[i] */
	    newton->F[i] = Fval;
	}
	    
	/* actual Jacobi loop to solve Jx = F */
	k = 0;
	while ((k < newton->jacobiIterations) && (resid > rtol)) {
	    //printf("Jacobi iteration %d, residual %.20f\n", k, r);
	    
	    /* xp = Dinv(b - Rx) */
	    for (i = 0; i < numjacobirows; i++) {
		row = newton->jacobiRows[i];
		val = newton->F[row];
		for (j = newton->J->rowStart[row]; j < newton->J->rowStart[row+1]; j++) {
		    val -= newton->J->values[j] * newton->x[newton->J->colIndex[j]];
		}
		val *= newton->d[row];
		newton->xp[row] = val;
	    }
	    
	    /* do a second iteration to get back from xp to x, and compute the relative residual */
	    resid = 0.0;
	    for (i = 0; i < numjacobirows; i++) {
		row = newton->jacobiRows[i];
		val = newton->F[row];
		for (j = newton->J->rowStart[row]; j < newton->J->rowStart[row+1]; j++) {
		    val -= newton->J->values[j] * newton->xp[newton->J->colIndex[j]];
		}
		val *= newton->d[row];
		newton->x[row] = val;
		
		/* compute relative residual */
		val -= newton->xp[row];
		resid += val*val;
	    }
	    
	    k += 2;
	}

	/* update r from solver result */
	val = 0.0;
	for (i = 0; i < len; i++) {
	    r[i] = r[i] - newton->x[i];
	    val += newton->x[i] * newton->x[i];
	}

	if (val < (newton->tolerance*newton->tolerance)) break;
    }

    /* update R: R = R + Q*r */
    for (i = 0; i < len; i++) {
	R[i] += q[i] * r[i];
    }
}



void interp1(double *x, double *v, double *xq, int lenx, int lenxq, double *result)
{
    int i, j;
    for (i = 0; i < lenxq; i++) {
	if (xq[i] < x[0]) {
	    // far LHS
	    result[i] = v[0];
	}
	else if (xq[i] >= x[lenx - 1]) {
	    // far RHS
	    result[i] = v[lenx - 1];
	}
	else {
	    // find the points it falls between
	    for (j = 0; j < (lenx-1); j++) {
		if ((x[j] <= xq[i]) && (x[j+1] > xq[i])) {
		    // interpolate from those points
		    double alpha = (xq[i] - x[j]) / (x[j+1] - x[j]);
		    result[i] = (1.0 - alpha) * v[j] + alpha * v[j+1];
		}
	    }
	}
    }
}

/*
 * As above, but x and v are interleaved in a single array
 */
void interp1_interleaved(double *xandv, double *xq, int lenx, int lenxq, double *result)
{
    int i, j;
    for (i = 0; i < lenxq; i++) {
	if (xq[i] < xandv[0]) {
	    // far LHS
	    result[i] = xandv[1];
	}
	else if (xq[i] >= xandv[(lenx - 1)*2]) {
	    // far RHS
	    result[i] = xandv[(lenx - 1)*2 + 1];
	}
	else {
	    // find the points it falls between
	    for (j = 0; j < (lenx-1); j++) {
		if ((xandv[j*2] <= xq[i]) && (xandv[(j+1)*2] > xq[i])) {
		    // interpolate from those points
		    double alpha = (xq[i] - xandv[j*2]) / (xandv[(j+1)*2] - xandv[j*2]);
		    result[i] = (1.0 - alpha) * xandv[j*2+1] + alpha * xandv[(j+1)*2+1];
		}
	    }
	}
    }
}

double factorial(double n)
{
    double i;
    double result = 1.0;
    for (i = 1.0; i <= n; i += 1.0) {
	result *= i;
    }
    return result;
}


bool croutDecomposition(double *A, double *L, double *U, int n)
{
    int i, j, k;
    double sum;

    for (i = 0; i < n; i++) {
	U[(i*n)+n] = 1.0;
    }

    for (j = 0; j < n; j++) {
	for (i = j; i < n; i++) {
	    sum = 0.0;
	    for (k = 0; k < j; k++) {
		sum = sum + L[(i*n)+k] * U[(k*n)+j];

	    }
	    L[(i*n)+j] = A[(i*n)+j] - sum;
	}
	for (i = j; i < n; i++) {
	    sum = 0.0;
	    for (k = 0; k < j; k++) {
		sum = sum + L[(j*n)+k] * U[(k*n)+i];
	    }
	    if (L[(j*n)+j] == 0.0) return false;
	    U[(j*n)+i] = (A[(j*n)+i] - sum) / L[(j*n)+j];
	}
    }
    return true;
}

void croutSolve(double *L, double *U, double *b, double *x, double *y, int n)
{
    int i, j;
    double sum;

    // forward solve
    for (i = 0; i < n; i++) {
	sum = 0.0;
	for (j = 0; j < i; j++) {
	    sum += L[(i*n)+j] * y[j];
	}
	y[i] = (b[i] - sum) / L[(i*n)+i];
    }

    // backward solve
    for (i = (n-1); i >= 0; i--) {
	sum = 0.0;
	for (j = (i+1); j < n; j++) {
	    sum += U[(i*n)+j] * x[j];
	}
	x[i] = (y[i] - sum); // U has unit diagonal so no need to divide
    }
}

/*
 * Gets the eigenvalues (in val) and the corresponding eigenvectors (in vec)
 * of a dense symmetric matrix A, size NxN, using the Jacobi method.
 * The eigenvalues are not sorted by this function.
 */
void getEigenvalues(int N, double *A, double *val, double *vec)
{
    int i, j, m, n;
    double threshold, theta, tau, sum;
    double p, q, r, s, t;
    double *tmp1, *tmp2;

    /* allocate temporary storage */
    tmp1 = new double[N];
    tmp2 = new double[N];

    /* initialise eigenvectors to identity matrix */
    for (m = 0; m < N; m++) {
	for (n = 0; n < N; n++) {
 	    if (m == n) vec[(m*N)+n] = 1.0;
	    else vec[(m*N)+n] = 0.0;
	}
    }
    for (m = 0; m < N; m++) {
	/* initialise to diagonal of matrix */
	tmp1[m] = A[(m*N)+m];
	val[m] = A[(m*N)+m];

	tmp2[m] = 0.0;
    }

    /* do Jacobi iterations */
    for (i = 0; i < 50; i++) {
	/* sum off-diagonal elements */
	sum = 0.0;
	for (m = 0; m < (N-1); m++) {
	    for (n = (m+1); n < N; n++) {
		sum += fabs(A[(m*N)+n]);
	    }
	}
	if (sum == 0.0) break; /* done! */

	/* set threshold */
	if (i < 5) threshold = 0.2*sum/((double)(N*N));
	else threshold = 0.0;

	for (m = 0; m < (N-1); m++) {
	    for (n = (m+1); n < N; n++) {
		if ((i >= 4) && (fabs(A[(m*N)+n]) == 0.0)) {
		    A[(m*N)+n] = 0.0;
		}
		else if (fabs(A[(m*N)+n]) > threshold) {
		    r = val[n] - val[m];
		    if (fabs(A[(m*N)+n]) == 0.0) {
			t = A[m*N+n] / r;
		    }
		    else {
			theta = 0.5*r/(A[m*N+n]);
			t = 1.0 / (fabs(theta) + sqrt(1.0+theta*theta));
			if (theta < 0.0) t = -t;
		    }
		    p = 1.0 / sqrt(1.0+t*t);
		    s = t*p;
		    tau = s / (1.0+p);
		    r = t * A[m*N+n];
		    tmp2[m] -= r;
		    val[m] -= r;
		    tmp2[n] += r;
		    val[n] += r;
		    A[m*N+n] = 0.0;

		    /* perform rotations */
		    for (j = 0; j <= (m-1); j++) {
			q=A[j*N+m];
			r=A[j*N+n];
			A[j*N+m]=q-s*(r+q*tau);
			A[j*N+n]=r+s*(q-r*tau);
		    }
		    for (j = m+1; j <= (n-1); j++) {
			q=A[m*N+j];
			r=A[j*N+n];
			A[m*N+j]=q-s*(r+q*tau);
			A[j*N+n]=r+s*(q-r*tau);
		    }
		    for (j = n+1; j < N; j++) {
			q=A[m*N+j];
			r=A[n*N+j];
			A[m*N+j]=q-s*(r+q*tau);
			A[n*N+j]=r+s*(q-r*tau);
		    }
		    for (j = 0; j < N; j++) {
			q=vec[j*N+m];
			r=vec[j*N+n];
			vec[j*N+m]=q-s*(r+q*tau);
			vec[j*N+n]=r+s*(q-r*tau);
		    }
		}
	    }
	}

	for (m = 0; m < N; m++) {
	    tmp1[m] += tmp2[m];
	    val[m] = tmp1[m];
	    tmp2[m] = 0.0;
	}
    }

    /* free temporary storage */
    delete[] tmp1;
    delete[] tmp2;
}

void denseMatrixVectorMultiply(double *A, double *x, double *b, int m, int n)
{
    int i, j;
    double val;

    /* loop over rows of matrix, and elements of result vector */
    for (i = 0; i < m; i++) {
	val = 0.0;
	/* loop over columns of matrix, and elements of source vector */
	for (j = 0; j < n; j++) {
	    val += x[j] * A[(j*m)+i];

	}
	b[i] = val;
    }
}

void denseMatrixVectorMultiplyTransposed(double *A, double *x, double *b, int m, int n)
{
    int i, j;
    double val;

    /* loop over rows of transposed matrix (columns of original), and result vector */
    for (i = 0; i < n; i++) {
	val = 0.0;
	/* loop over columns of transposed matrix (rows of original), and source vector */
	for (j = 0; j < m; j++) {
	    val += x[j] * A[(i*m)+j];
	}
	b[i] = val;
    }
}

void denseCholeskyDecomp(int N, double *B, double *L, double *U)
{
    int i, j, k;
    double sum;

    memcpy(L, B, N*N*sizeof(double));

    for (i = 0; i < N; i++) {
	for (j = i; j < N; j++) {
	    sum = L[(i*N)+j];
	    for (k = (i-1); k >= 0; k--) {
		sum -= L[(i*N)+k] * L[(j*N)+k];
	    }
	    if (i == j) {
		if (sum <= 0.0) {
		    // matrix is not positive definite
		    return;
		}
		L[(i*N)+i] = sqrt(sum);
		U[(i*N)+i] = sqrt(sum);
	    }
	    else {
		L[(j*N)+i] = sum / L[(i*N)+i];
		U[(i*N)+j] = sum / L[(i*N)+i];
		
		/* make sure other elements are zeroed */
		L[(i*N)+j] = 0.0;
		U[(j*N)+i] = 0.0;
	    }
	}
    }
}

void invertLowerTriangle(int N, double *L, double *I)
{
    int i, j, k;
    double sum;

    memcpy(I, L, N*N*sizeof(double));

    for (i = 0; i < N; i++) {
	I[(i*N)+i] = 1.0 / L[(i*N)+i];
	for (j = i+1; j < N; j++) {
	    sum = 0.0;
	    for (k = i; k < j; k++) {
		sum -= I[(j*N)+k] * I[(k*N)+i];
	    }
	    I[(j*N)+i] = sum / L[(j*N)+j];
	}
    }
}

void transposeDenseMatrix(int N, double *in, double *out)
{
    int i, j;
    for (i = 0; i < N; i++) {
	for (j = 0; j < N; j++) {
	    out[(j*N)+i] = in[(i*N)+j];
	}
    }
}

void transposeMatrixRectangular(int Nr, int Nc, double *in, double *out)
{
    int i, j;
    for (i = 0; i < Nc; i++) {
	for (j = 0; j < Nr; j++) {
	    out[(j*Nc)+i] = in[(i*Nr)+j];
	}
    }
}

void denseMatrixMatrixMultiply(int N, double *A, double *B, double *out)
{
    /* out(i,j) = row i in A * column j in B */
    int i, j, k;
    double sum;
    for (i = 0; i < N; i++) {
	for (j = 0; j < N; j++) {
	    sum = 0.0;
	    for (k = 0; k < N; k++) {
		sum += A[(i*N)+k] * B[(k*N)+j];
	    }
	    out[(i*N)+j] = sum;
	}
    }
}

void denseMatrixMatrixMultiplyRectangular(int Nr, int Nc, double *A, double *B, double *out)
{
    /* out(i,j) = row i in A * column j in B */
    int i, j, k;
    double sum;
    for (i = 0; i < Nr; i++) {
	for (j = 0; j < Nr; j++) {
	    sum = 0.0;
	    for (k = 0; k < Nc; k++) {
		sum += A[(i*Nc)+k] * B[(k*Nr)+j];
	    }
	    out[(i*Nr)+j] = sum;
	}
    }
}



double Trapz2Dpolar(double* f , double r_min, double r_max, double th_min, double th_max, int Nr, int Nth ) {

	double I=0;
	double I2 = 0;
	double I4 = 0;
	double r;

	double dr, dth;

	int i,j;
	

	dr = (r_max-r_min)/(Nr-1);

	dth = (th_max-th_min)/(Nth-1);
	
	//I = r(1)*(f(1,1)+f(1,end))+r(end)*(f(end,1)+f(end,end));
	I = r_min*(f[0]+f[Nth-1])+r_max*(f[(Nr-1)*Nth]+f[(Nr-1)*Nth+(Nth-1)]);
	

	for (i = 1; i<(Nr-1); i++){
		r = r_min+dr*i;
	    I2 = I2 + r*(f[i*Nth] + f[i*Nth+(Nth-1)]) ;

	    for (j = 1; j<(Nth-1); j++){
	        I4 = I4 + r*f[i*Nth+j] ;
	    }

	}

	for (j = 1; j<(Nth-1); j++){
	    I2 = I2 + r_min*f[j] + r_max*f[(Nr-1)*Nth+j];
	}

	I = (I + 2*I2 + 4*I4)*0.25*dr*dth;

	return I;
}

double Trapz1DPolar(double* f , double r_min, double r_max, int Nr ) {

	double I=0;
	double I2 = 0;
	double r;

	double dr;

	

	dr = (r_max-r_min)/(Nr-1);

	I = r_min*f[0] + r_max*f[(Nr-1)];
	

	for (int i = 1; i<(Nr-1); i++){
		r = r_min+dr*i;
	    I2 = I2 + r*f[i];
	}

	
	I = (I + 2*I2)*0.5*dr;

	return I;
}

double ScalarProductPolar(double *f1, double *f2, int Nr, int Nth){


	int i;
	double sp;

	double *f = new double [Nr*Nth];

	for (i=0; i<(Nr*Nth);i++){

			f[i] = f1[i]*f2[i];		

	}
	

	sp = Trapz2Dpolar(f , 0, 1, 0, 2*M_PI, Nr, Nth);

	delete[] f;

	return sp;

}

double Trapz2DCartesian(double* f , double x_min, double x_max, double y_min, double y_max, int Nx, int Ny ) 
{

	double I=0;
	double I2 = 0;
	double I4 = 0;


	double dx, dy;

	int i,j;
	

	dx = (x_max-x_min)/(Nx-1);

	dy = (y_max-y_min)/(Ny-1);

	
	//I = r(1)*(f(1,1)+f(1,end))+r(end)*(f(end,1)+f(end,end));
	I = (f[0]+f[Ny-1]) + (f[(Nx-1)*Ny]+f[(Nx-1)*Ny+(Ny-1)]);
	

	for (i = 1; i<(Nx-1); i++){
		
	    I2 = I2 + (f[i*Ny] + f[i*Ny+(Ny-1)]) ;

	    for (j = 1; j<(Ny-1); j++){
	        I4 = I4 + f[i*Ny+j] ;
	    }

	}

	for (j = 1; j<(Ny-1); j++){
	    I2 = I2 + f[j] + f[(Nx-1)*Ny+j];
	}

	I = (I + 2*I2 + 4*I4)*0.25*dx*dy;

	return I;
}


double Trapz1DCartesian(double* f , double x_min, double x_max, int pointsx ) 
{

	double I=0;
	double I2 = 0;
	double x;

	double dx;

	

	dx = (x_max-x_min)/(pointsx-1);

	I = f[0] + f[(pointsx-1)];
	

	for (int i = 1; i<(pointsx-1); i++){
		x = x_min+dx*i;
	    I2 = I2 + f[i];
	}

	
	I = (I + 2*I2)*0.5*dx;

	return I;
}

double ScalarProductCartesian(double *f1, double *f2, double Lx, double Ly, int Nx, int Ny)
{

	int i;
	double sp;

	double *f = new double [Nx*Ny];

	for (i=0; i<(Nx*Ny);i++){

			f[i] = f1[i]*f2[i];

	}


	sp = Trapz2DCartesian(f , 0, Lx, 0, Ly, Nx, Ny);

	delete[] f;

	return sp;

}

//bool gaussSolver(double* mat, double* rhs, double* minPivot, int* row, int nbRows)
bool gaussSolver(double* mat, double* rhs, int nbRows)
{
  // returns false when solver fails; use last args to display error
  //
  // minPivot is absolute magnitude of minimal pivot
  // row is row on which a null pivot is found when algorithm fails
  // (in other words row is the matrix rank)
  // is determinant value required on exit ?

  bool flag = false; // useless unless rhs.size() = 0 !!
  int  p_row=0, p_r=0, row_i=0;
  double det=1.0, sw_m=0.0, pivot=0.0, mult=0.0;
  double abs_pivot=0.0;
  double theZeroThreshold = 1e-100;
  int row = 0;
  double minPivot = 1.e+30;
  int rc = 0, p_rc = 0, rr = 0, p_rr = 0;

  
  //
  // Gaussian elimination (matrix triangulation)
  //
  for (row = 0; row < nbRows; row++)
    {
      // current step number
      p_row = row;
      p_rc = nbRows * row + row;
      pivot = mat[p_rc];
      abs_pivot = std::abs(pivot);
      // search for partial pivot of maximum absolute magnitude in row-th column
      rc = p_rc + nbRows;
      for (row_i = row + 1; row_i < nbRows; rc += nbRows, row_i++)
        {
          if (abs_pivot < std::abs(mat[rc]))
            {
              p_row = row_i;
              pivot = mat[rc];
              abs_pivot = std::abs(pivot);
            }
        }
      if (abs_pivot < minPivot)
        {
          minPivot = abs_pivot;
          if (abs_pivot < theZeroThreshold) return false; // non invertible matrix
        }
      // swap entries on current row and pivot row as well as on right hand side (if needed)
      if (p_row > row)
        {
          rc = nbRows * row + row;
          p_rc = nbRows * p_row + row;

          for (int col = row; col < nbRows; rc++, p_rc++, col++)
            {
              sw_m = mat[rc];
              mat[rc] = mat[p_rc];
              mat[p_rc] = sw_m;
            }
          sw_m = rhs[row];
          rhs[row] = rhs[p_row];
          rhs[p_row] = sw_m;
        }

      // Gaussian elimination in rows row+1, ..., nbRows
      det *= pivot;
      pivot = 1. / pivot;
      p_rr = nbRows * row + row;
      rr = p_rr + nbRows;
      for (row_i = row + 1; row_i < nbRows; row_i++, rr += nbRows)
        {
          p_rc = p_rr;
          rc = rr;
          mult = - pivot * mat[rc];
          for (int col = row; col < nbRows; p_rc++, rc++, col++)
            mat[rc] += mult * mat[p_rc];
          rhs[row_i] += mult * rhs[row];
        }
      flag = true;
    }

  // solve upper triangularized matrix : backward substitution
  for (p_r = nbRows; p_r > 0; p_r--)
    {
      row = p_r - 1;
      rc = nbRows * row + row;
      p_rc = rc + 1;
      for (int col = row + 1; col < nbRows; col++, p_rc++)
        rhs[row] -= mat[p_rc] * rhs[col];
      rhs[row] /= mat[rc];
    }
  return flag;
}


