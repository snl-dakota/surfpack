/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#ifndef SURFPACK_SYSTEM_HEADERS_H
#define SURFPACK_SYSTEM_HEADERS_H

#ifdef HAVE_CONFIG_H
#include "surfpack_config.h"
#else
#include "surf77_config.h"
#endif

#include <algorithm>
#include <limits>
#include <cassert>
#include <cctype>
#include <cerrno>
#include <cfloat>
#include <climits>
#include <cmath>
#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iomanip> 
#include <iostream>
#include <iterator>
#include <list>
#include <map>
#include <numeric>
#include <queue>
#include <set>
#include <sstream>
#include <stack>
#include <stdexcept>
#include <string>
#include <sys/time.h>
#include <vector>

typedef std::vector<double> VecDbl;
typedef std::vector<double>::const_iterator VecDblIt;
typedef std::vector<unsigned> VecUns;
typedef std::vector<unsigned>::const_iterator VecUnsIt;
typedef std::vector< std::vector< unsigned > > VecVecUns;
typedef std::vector< std::vector< double > > VecVecDbl;
typedef std::vector< std::string > VecStr;
typedef std::set< unsigned > SetUns;

typedef std::pair< std::string, std::string > ModelParam;
typedef std::map< std::string, std::string> ParamMap;
typedef std::pair< std::string, ParamMap> Command;


// Commonly used BLAS/LAPACK subroutines

#ifdef HAVE_CONFIG_H
// Tolerate F77_FUNC macro redefinition warnings in the autotools build
#define DGETRF_F77 F77_FUNC(dgetrf,DGETRF)
#define DGETRI_F77 F77_FUNC(dgetri,DGETRI)
#define DGEMV_F77  F77_FUNC(dgemv,DGEMV)
#define DGEMM_F77  F77_FUNC(dgemm,DGEMM)
#define DDOT_F77   F77_FUNC(ddot, DDOT)
#define DGELS_F77  F77_FUNC(dgels,DGELS)
#define DGGLSE_F77 F77_FUNC(dgglse,DGGLSE)

#else
// Use the CMake generated fortran name mangling macros (eliminate warnings)
#define DGETRF_F77 SURF77_GLOBAL(dgetrf,DGETRF)
#define DGETRI_F77 SURF77_GLOBAL(dgetri,DGETRI)
#define DGEMV_F77  SURF77_GLOBAL(dgemv,DGEMV)
#define DGEMM_F77  SURF77_GLOBAL(dgemm,DGEMM)
#define DDOT_F77   SURF77_GLOBAL(ddot, DDOT)
#define DGELS_F77  SURF77_GLOBAL(dgels,DGELS)
#define DGGLSE_F77 SURF77_GLOBAL(dgglse,DGGLSE)
#endif

#ifdef __cplusplus
extern "C" { // prevent C++ name mangling
#endif

// Perform LU factorization
void DGETRF_F77(const int* m, const int* n, double* a, const int* lda,
                int* ipiv, int* info);

// Compute the inverse of a matrix expressed as an LU decomposition
// (i.e., call dgetrf on the matrix first)
void DGETRI_F77(const int* n, double* a, const int* lda, const int* ipiv,
                double* work, const int* lwork, int* info);

// Matrix-vector multiplication
void DGEMV_F77(char* trans, const int* m, const int* n, const double* alpha,
               const double* A, const int* lda, const double* x,
               const int* incx, const double* beta, double* y, const int* incy);

// Matrix-matrix multiplication
void DGEMM_F77(char* transa, char* transb, const int* m, const int* n,
               const int* k, const double* alpha, const double* A,
               const int* lda, const double* B, const int* ldb,
               const double* beta, double* C, const int* ldc);

// Vector-vector inner product
double DDOT_F77(const int* n, const double* x, const int* incx,
                const double* y, const int* incy);

// Least-squares solution to linear system of equations
void DGELS_F77(const char* trans, const int* nrows, const int* ncols,
               const int* nrhs, double* A, const int* lda, double* b,
               const int* ldb, double* work, const int* lwork, int* info);

// Performs least-squares solve subject to equality constraints
void DGGLSE_F77(const int* m, const int* n, const int* p, double* A,
                const int* lda, double* B, const int* ldb, double* c,
                double* d, double* x, double* work, const int* lwork,
                int* info);

}  // end extern "C" of BLAS/LAPACK subroutines

#endif // SURFPACK_SYSTEM_HEADERS_H

