/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#ifndef __SURFPACK_LAPACK_WRAPPERS_H__
#define __SURFPACK_LAPACK_WRAPPERS_H__

#ifdef HAVE_CONFIG_H
#include "surfpack_config.h"
/* WJB - ToDo:  one more iteration to get the CMake build functional
#elif HAVE_EMPTY_CONFIG_H
#include "surf77_config.h"
*/
#endif

/***************************************************************************/
/***************************************************************************/
/**** The BLAS and LAPACK wrappers should be the same whichever version ****/
/**** of SurfMat we use, and since they are _ONLY_ wrappers they should ****/
/**** the should be inline functions so they should be in a header file ****/
/**** so it makes sense to put them here.                               ****/
/***************************************************************************/
/***************************************************************************/

#ifdef HAVE_CONFIG_H
#define DGETRF_F77 F77_FUNC(dgetrf,DGETRF)
#define DGETRI_F77 F77_FUNC(dgetri,DGETRI)
#define DGEMV_F77  F77_FUNC(dgemv,DGEMV)
#define DGEMM_F77  F77_FUNC(dgemm,DGEMM)
#define DDOT_F77   F77_FUNC(ddot, DDOT)
#define DGELS_F77  F77_FUNC(dgels,DGELS)

#define DPOTRF_F77 F77_FUNC(dpotrf,DPOTRF)
#define DPOTRI_F77 F77_FUNC(dpotri,DPOTRI)
#define DPOTRS_F77 F77_FUNC(dpotrs,DPOTRS)
#define DPOCON_F77 F77_FUNC(dpocon,DPOCON)
#define DGETRS_F77 F77_FUNC(dgetrs,DGETRS)
#define DLANGE_F77 F77_FUNC(dlange,DLANGE)
#define DGECON_F77 F77_FUNC(dgecon,DGECON)
#define DGGLSE_F77 F77_FUNC(dgglse,DGGLSE)
#define DSYEV_F77  F77_FUNC(dsyev,DSYEV)
/*
#elif HAVE_EMPTY_CONFIG_H
// Use the CMake generated fortran name mangling macros (eliminate warnings)
#define DGETRF_F77 SURF77_GLOBAL(dgetrf,DGETRF) 
#define DGETRI_F77 SURF77_GLOBAL(dgetri,DGETRI) 
#define DGEMV_F77  SURF77_GLOBAL(dgemv,DGEMV) 
#define DGEMM_F77  SURF77_GLOBAL(dgemm,DGEMM) 
#define DDOT_F77   SURF77_GLOBAL(ddot, DDOT) 
#define DGELS_F77  SURF77_GLOBAL(dgels,DGELS)

#define DPOTRF_F77 SURF77_GLOBAL(dpotrf,DPOTRF)
#define DPOTRI_F77 SURF77_GLOBAL(dpotri,DPOTRI)
#define DPOTRS_F77 SURF77_GLOBAL(dpotrs,DPOTRS)
#define DPOCON_F77 SURF77_GLOBAL(dpocon,DPOCON)
#define DGETRS_F77 SURF77_GLOBAL(dgetrs,DGETRS)
#define DLANGE_F77 SURF77_GLOBAL(dlange,DLANGE)
#define DGECON_F77 SURF77_GLOBAL(dgecon,DGECON)
#define DGGLSE_F77 SURF77_GLOBAL(dgglse,DGGLSE)
#define DSYEV_F77  SURF77_GLOBAL(dsyev,DSYEV)
*/
#endif

/***************************************************************************/
/**** BLAS Fortran to C name mangling                                   ****/
/***************************************************************************/

extern "C" { // prevent C++ name mangling

// Vector-vector inner product
double DDOT_F77(const int* n, const double* x, const int* incx,
		const double* y, const int* incy);


// Matrix-vector multiplication
void DGEMV_F77(char* trans, const int* m, const int* n, const double* alpha, 
	       const double* A, const int* lda, const double* x,
	       const int* incx, const double* beta, double* y, const int* incy);

// Matrix-matrix multiplication
void DGEMM_F77(char* transa, char* transb, const int* m, const int* n,
	       const int* k, const double* alpha, const double* A,
	       const int* lda, const double* B, const int* ldb, 
	       const double* beta, double* C, const int* ldc);

/***************************************************************************/
/**** LAPACK Fortran to C name mangling                                 ****/
/***************************************************************************/

// Perform Cholesky factorization
void DPOTRF_F77(const char* uplo, const int* n, double* AChol, const int* lda, int* info);


// Compute the inverse of a matrix expressed as an cholesky decomposition (i.e., call dpotrf on the matrix first)
void DPOTRI_F77(const char* uplo, const int* n, double* ACholInv, const int* lda, int* info);

// solve A*X=B for X, where A={A || A^T} after A has been Cholesky factorized (i.e., call dptorf on the matrix first)
void DPOTRS_F77(const char* uplo, const int* n, const int* nRHS, 
		const double* AChol,
		const int* ldAChol , double* RHS, 
		const int* ldRHS, int* info);


//function to compute the condition number of a matrix from the Cholesky factorization
void DPOCON_F77(const char* uplo, const int* n, const double* AChol, const int* lda, const double* anorm, double* rconda, double* work, int* iwork, int* info);
  
  




// Perform LU factorization
void DGETRF_F77(const int* m, const int* n, double* a, const int* lda,
		int* ipiv, int* info);


// Compute the inverse of a matrix expressed as an LU decomposition (i.e., call dgetrf on the matrix first)
void DGETRI_F77(const int* n, double* a, const int* lda, const int* ipiv,
		double* work, const int* lwork, int* info);


// solve A*X=B for X, where A={A || A^T} after A has been LU factorize (i.e., call dgetrf on the matrix first)
void DGETRS_F77(char* transLU, const int* n, const int* nRHS, 
		const double* LU, const int* ldLU , const int* ipiv, 
		double* RHS, const int* ldRHS, int* info);

//function to compute the norm of a matrix A, choices are
//M max(abs(A(i,j))) this is not a consistent matrix norm, 
//1 one norm of a matrix, maximum column sum, 
//I infinity norm of matrix, maximum row sum,or 
//F frobenius norm of a matrix, square root of sum of squares
double DLANGE_F77(char *whichnorm, int *M, int *N, const double *A, int *LDA,
		  double *work);

//function to compute the condition number of a matrix
void DGECON_F77(char *whichnorm, int *N, const double *ALU, int *LDA, double *anorm,
		double *rcond, double *work, int *iwork, int *info);

// Least-squares solution to linear system of equations
void DGELS_F77(const char* trans, const int* nrows, const int* ncols,
	       const int* nrhs, double* A, const int* lda, double* b,
	       const int* ldb, double* work, const int* lwork, int* info);

// Performs least-squares solve subject to equality constraints
void DGGLSE_F77(const int* m, const int* n, const int* p, double* A,
		const int* lda, double* B, const int* ldb, double* c,
		double* d, double* x, double* work, const int* lwork,
		int* info);

// determines eigenvalues and (optionally) eigenvectors for a real symmetric matrix
void DSYEV_F77(const char* jobz, const char* uplo, const int* N, 
	       double *A_EIGVECT, const int* lda, double* eigval, 
	       double* work, const int* lwork, int* info);

} // extern "C" (prevent C++ name mangling)


#endif // __SURFPACK_LAPACK_WRAPPERS_H__

