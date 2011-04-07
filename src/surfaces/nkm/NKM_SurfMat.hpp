#ifndef __SURFMAT_HPP__
//do not #define __SURFMAT_HPP__ here, that should/must only be done in either CustomSurfMat.hpp OR TeuchosSurfMat.hpp, to keep them from being included directly when the other has already been included

//#define __SURFMAT_ERR_CHECK__
//#define __SURFMAT_ZERO_MEM__

#ifdef HAVE_CONFIG_H
#include "surfpack_config.h"
#endif

#include "surfpack_LAPACK_wrappers.h"

#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cmath>


#include "NKM_SurfMat_Native.hpp" //a native implementation
//#include "SurfMat_Teuchos.hpp" //a wrapper for the Teuchos Serial Dense Matrix class


namespace nkm {

typedef SurfMat<double> MtxDbl;
typedef SurfMat<int> MtxInt;


/***************************************************************************/
/**** The BLAS wrappers start here                                      ****/ 
/***************************************************************************/

/// sum of element by element products of a and b (works for matrices as well as vectors and that functionality is required for the Kriging/GP implementation), wraps DDOT
inline double dot_product(const MtxDbl& a, const MtxDbl& b)
{
  int nelem=a.getNRows()*a.getNCols();
#ifdef __SURFMAT_ERR_CHECK__
  assert(nelem==b.getNRows()*b.getNCols());
#endif
  int inc=1;
  // ddot will not violate the constness
  return DDOT_F77(&nelem, a.ptr(0), &inc, b.ptr(0), &inc);
}


/// matrix matrix Mult OR Matrix Vector Mult: C=scaleAB*A*B+scaleC*C, where A={A || A^T} and B={B || B^T}, wraps and chooses between DGEMV and DGEMM
inline MtxDbl& matrix_mult(MtxDbl& C, const MtxDbl& A, const MtxDbl& B, 
			   double scaleC=0.0, double scaleAB=1.0, 
			   char transA='N', char transB='N')
{
#ifdef __SURFMAT_ERR_CHECK__
  assert(((transA=='N')||(transA=='T'))&&((transB=='N')||(transB=='T')));  
#endif
  int nrowsC, ncolsC, ninnerA, ninnerB;
  if(transA=='N') {
    nrowsC =A.getNRows();
    ninnerA=A.getNCols();
  }
  else{
    nrowsC =A.getNCols();
    ninnerA=A.getNRows();
  }
  if(transB=='N') {
    ncolsC =B.getNCols();
    ninnerB=B.getNRows();
  }
  else{
    ncolsC =B.getNRows();
    ninnerB=B.getNCols();
  }
#ifdef __SURFMAT_ERR_CHECK__
  if(!(ninnerA==ninnerB)){
    printf("ninnerA=%d ninnerB=%d",ninnerA,ninnerB);
    printf("\n");
    assert(ninnerA==ninnerB);
  }
#endif
  C.newSize(nrowsC,ncolsC);
  C.putTol(A.getTol()); //revise this once nonzero-tol starts to be used
  
  int nrowsA= A.getNRows(); // leading dimension irrespective of trans
  int ncolsA= A.getNCols();
  int nrowsB= B.getNRows(); // leading dimension irrespective of trans
  int inc=1;
  if(ncolsC==1)//BLAS2 matrix vector multiply
    { 
      
      //printf("transA=%c transB=%c nrowsC=%d ninnnerA=%d ninnerB=%d A.getNRows()=%d\n",transA,transB,nrowsC,ninnerA,ninnerB,A.getNRows());
      
      //DGEMV_F77(&transA,&nrowsC,&ninnerA,&scaleAB,A.ptr(0,0),&nrowsC,B.ptr(0),&inc,&scaleC,C.ptr(0),&inc); 
      DGEMV_F77(&transA,&nrowsA,&ncolsA,&scaleAB,A.ptr(0,0),&nrowsA,B.ptr(0),&inc,&scaleC,C.ptr(0),&inc);
      
    }
  else //BLAS3 matrix matrix multiply
    DGEMM_F77(&transA,&transB,&nrowsC,&ncolsC,&ninnerA,&scaleAB,A.ptr(0,0),&nrowsA,B.ptr(0,0),
	      &nrowsB,&scaleC,C.ptr(0,0),&nrowsC);
  return C;
}


/***************************************************************************/
/**** The LAPACK wrappers start here                                    ****/ 
/***************************************************************************/

/// perform a numerically optimally preconditioned Cholesky factorization of a real symmetric positive semi-definite matrix, wraps DPOTRF, returns the lower triangular portion, the preconditioning is "undone" so the L matrix is the same (minus some rounding error due to poor conditioning) as it would be without precondtioning, if info>0 then the preconditined matrix is singular, rcondprecond is the standard fast estimate of the reciprocal of the condition number of the numerically optimally preconditioned real symmetric positive definite matrix (that's what affects the round off error)
inline MtxDbl& Chol_fact(MtxDbl& matrix, int& info, double& rcondprecond)
{
  int nrows = static_cast<int>(matrix.getNRows());
  int ncols = static_cast<int>(matrix.getNCols());
#ifdef __SURFMAT_ERR_CHECK__
  assert(nrows==ncols);
#endif
  //printf("Cholf_fact() nrows=%d ncols=%d\n",nrows,ncols);
  char uplo='L';
  int lda=nrows;
  int info_local=0;

  MtxDbl work(3*nrows);
  MtxInt iwork(nrows);

  int abspower;
  int minabspower;
  int maxabspower;
  MtxDbl scalefactor(nrows);
  double log_of_2=std::log(2.0);
  abspower=std::floor(0.5+std::log(std::sqrt(matrix(0,0)))/log_of_2);
  scalefactor(0)=std::pow(2.0,(double) -abspower);
  //abspower=log2(sqrt(matrix(0,0)));
  //scalefactor(0)=1.0/sqrt(matrix(0,0));
  minabspower=maxabspower=abspower;
  for(int i=1; i<nrows; ++i) {
    abspower=std::floor(0.5+std::log(std::sqrt(matrix(i,i)))/log_of_2);
    scalefactor(i)=std::pow(2.0,(double) -abspower); //this is the "numerically optimal" preconditioning of a real symmetric positive definite matrix by "numerically optimal" I meant the analytically optimal (for reducing condition number) scaling has been rounded to the nearest power of 2 so that we don't lose any bits of accuracy due to rounding error due to preconditioning
    //abspower=log2(sqrt(matrix(i,i)));
    //scalefactor(i)=1.0/sqrt(matrix(i,i));
    minabspower=(abspower<minabspower)?abspower:minabspower;
    maxabspower=(abspower>maxabspower)?abspower:maxabspower;   
  }
  
  if(maxabspower!=minabspower) {
    //only do the preconditioning if the maximum and minimum numerically optimal scaling factors are different (by a factor of 2 or more) because otherwise the real symmetric positive definite matrix is already numerically optimally preconditioned.
    for(int j=0; j<nrows; ++j)
      for(int i=0; i<nrows; ++i)
	matrix(i,j)*=scalefactor(i)*scalefactor(j);
  }

  //calculate orignorm now so we don't have to store a copy of the unfactored matrix
  char whichnorm='1';
  double orignorm=DLANGE_F77(&whichnorm,&nrows,&ncols,matrix.ptr(0),&lda,work.ptr(0));

  DPOTRF_F77(&uplo,&nrows,matrix.ptr(0),&lda,&info_local);
#ifdef __SURFMAT_ERR_CHECK__
  assert(info_local>=0);
#endif
  info=info_local;

  //the rcond of the "numerically optimally" preconditioned real symmetric positive definte matrix, feed it orignorm
  DPOCON_F77(&uplo,&nrows,matrix.ptr(0),&lda,&orignorm,&rcondprecond,work.ptr(0),iwork.ptr(0),&info_local);
#ifdef __SURFMAT_ERR_CHECK__
  assert(info_local==0);
#endif

  if(maxabspower!=minabspower) {
    //undo the "numerically optimal" preconditioning of the real symmetric positive definite matrix, that is other than possibly avoiding rounding error due to poorly condition matrix this function produces the same cholesky "L" matrix that you would get without preconditiong.
    for(int i=0; i<nrows; ++i)
      scalefactor(i)=1.0/scalefactor(i); //multiplication can be faster than division

    for(int j=0; j<nrows; ++j)
      for(int i=j; i<nrows; ++i)  //it's lower triangular
	matrix(i,j)*=scalefactor(i);
  }

  return matrix;
}

/// inverts a real symmetric positive definite matrix in place after the Cholesky factorization has already been done, requires the lower triangular portion as input returns the full symmetric inverse (as opposed to just the lower triangular part), wraps DPOTRF
inline MtxDbl& inverse_after_Chol_fact(MtxDbl& matrix)
{
  int nrows = static_cast<int>(matrix.getNRows());
  int ncols = static_cast<int>(matrix.getNCols());
#ifdef __SURFMAT_ERR_CHECK__
  assert(nrows==ncols);
#endif
  char uplo='L';
  int lda=nrows;
  int info=0;

  DPOTRI_F77(&uplo,&nrows,matrix.ptr(0),&lda,&info);
#ifdef __SURFMAT_ERR_CHECK__
  assert(info==0);
#endif
  //fill in the top half of the inverse
  for(int j=0; j<ncols-1; ++j)
    for(int i=j+1; i<nrows; ++i)
      matrix(j,i)=matrix(i,j);

  return matrix;
}

//computes the (1 norm) reciprocal of the condition number of A from the Cholesky factorization of A (A must be real symmetric and positive semi-definite)
inline double rcond_after_Chol_fact(const MtxDbl& A, const MtxDbl& AChol)
{
  double rconda;
  char whichnorm='1';
  char uplo='L';
  int nrows=A.getNRows();
  int ncols=A.getNCols();
  int lda=AChol.getNRows();
  MtxDbl work(3*nrows);
  MtxInt iwork(nrows);
  int info;
#ifdef __SURFMAT_ERR_CHECK__
  assert((lda==nrows)&&(lda==AChol.getNCols())&&(nrows==ncols));
#endif
  double anorm=DLANGE_F77(&whichnorm,&nrows,&ncols,A.ptr(0),&lda,work.ptr(0));
  DPOCON_F77(&uplo,&nrows,AChol.ptr(0),&lda,&anorm,&rconda,work.ptr(0),iwork.ptr(0),&info);
#ifdef __SURFMAT_ERR_CHECK__
  assert(info==0);
#endif
  return rconda;
}

/// solves A*X=B for X, where A is symmetric positive definite and B={B || B^T}, without changing the contents of B, after A has been Cholesky factorized, AChol must contain the lower triangular portion of the factorization of A, wraps DPOTRS 
inline MtxDbl& solve_after_Chol_fact(MtxDbl& result, const MtxDbl& AChol, const MtxDbl& BRHS,char transB='N')
{
  int n   = static_cast<int>(AChol.getNRows());
#ifdef __SURFMAT_ERR_CHECK__
  assert(((transB=='N') && (BRHS.getNRows()==n)) ||
	 ((transB=='T') && (BRHS.getNCols()==n)));
#endif
  int lda = n;
  int ldb = n;
  char uplo='L';
  if(transB=='N'){ //copy B into the work matrix "result"
    /*
    result.newSize(BRHS.getNRows(),BRHS.getNCols());
    result.putTol(BRHS.getTol());
    for(int j=0; j<BRHS.getNCols(); j++)
      for(int i=0; i<BRHS.getNRows(); i++)
	result(i,j)=BRHS(i,j);
    */
    result = BRHS;
  }
  else{ //copy the transpose of B into work matrix "result"
    result.newSize(BRHS.getNCols(),BRHS.getNRows());
    result.putTol(BRHS.getTol());
    for(int i=0; i<BRHS.getNRows(); i++)
      for(int j=0; j<BRHS.getNCols(); j++)
	result(j,i)=BRHS(i,j);    
  }
  int nrhs= static_cast<int>(result.getNCols());
  
  //printf("solve_after_Chol_fact: n=%d nrhs=%d\n",n,nrhs);
  int info=0;
  DPOTRS_F77(&uplo, &n, &nrhs, AChol.ptr(0,0), &lda, result.ptr(0,0), &ldb, &info);
  return result;
}

/// perform an LU factorization with partial pivoting, wraps DGETRF
inline MtxDbl& LU_fact(MtxDbl& matrix, MtxInt& ipvt)
{
  int nrows = static_cast<int>(matrix.getNRows());
  int ncols = static_cast<int>(matrix.getNCols());
#ifdef __SURFMAT_ERR_CHECK__
  assert(nrows==ncols);
#endif
  // dgetrf may reorder the rows, the mapping between old and new rows
  // is returned in ipvt.
  if((ipvt.getNRows()!=nrows)||(ipvt.getNCols()!=1))
    ipvt.newSize(nrows); 
  
  int lda = ncols;
  int info = 0;
  //printf("nrows=%d ncols=%d\n",nrows,ncols);
  //std::cout << "Matrix size: " << n_rows << " " << n_cols << std::endl;
  DGETRF_F77(&nrows,&ncols,matrix.ptr(0,0),&lda,ipvt.ptr(0),&info);
  //std::cout << "Done with dgetrf" << std::endl;
  return matrix;
}

/// inverts a matrix, in place, after an LU factorization has already been done, wraps DGETRI
inline MtxDbl& inverse_after_LU_fact(MtxDbl& matrix, MtxInt& ipvt)
{
  int nrows = static_cast<int>(matrix.getNRows());
  int ncols = static_cast<int>(matrix.getNCols());
#ifdef __SURFMAT_ERR_CHECK__
  assert(nrows==ncols);
#endif
  int lwork = ncols;  // should be optimal blocksize
  MtxDbl work(lwork);
  int lda = nrows;
  int info = 0;
  //std::cout << "Matrix size: " << n_rows << " " << n_cols << std::endl;
  DGETRI_F77(&nrows,matrix.ptr(0,0),&lda,ipvt.ptr(0),work.ptr(0),&lwork,&info);
  //std::cout << "Done with getri" << std::endl;
  return matrix;
}



//computes the (1 norm) reciprocal of the condition number of A from the LU factorization of A
inline double rcond_after_LU_fact(const MtxDbl& A, const MtxDbl& ALU) {
  double rcond;
  char whichnorm='1';
  int M=A.getNRows();
  int LDA=ALU.getNRows();
  int N=ALU.getNCols();
  MtxDbl work(4*N); 
  MtxInt iwork(N);
  int info;
#ifdef __SURFMAT_ERR_CHECK__
  assert((LDA==N)&&(LDA==M)&&(LDA==A.getNCols()));
#endif
  double anorm=DLANGE_F77(&whichnorm,&M,&N,A.ptr(0,0),&LDA,work.ptr(0));
  DGECON_F77(&whichnorm,&N,ALU.ptr(0,0),&LDA,&anorm,&rcond,work.ptr(0),
	     iwork.ptr(0),&info);
  return rcond;
}

/// solves A*X=B for X, where A={A || A^T} and B={B || B^T}, without changing the contents of B, after A has been LU factorized, wraps DGETRS 
inline MtxDbl& solve_after_LU_fact(MtxDbl& result, const MtxDbl& ALU, const MtxInt& ipvt, 
				   const MtxDbl& BRHS, char transA='N', 
				   char transB='N')
{
  int n   = static_cast<int>(ALU.getNRows());
#ifdef __SURFMAT_ERR_CHECK__
  assert((transA=='N' || transA=='T')&&(n == ALU.getNCols())&&
	 (n=ipvt.getNRows()));
  assert(((transB=='N') && (BRHS.getNRows()==n)) ||
	 ((transB=='T') && (BRHS.getNCols()==n)));
#endif
  int lda = n;
  int ldb = n;
  if(transB=='N'){ //copy B into the work matrix "result"
    /*
    result.newSize(BRHS.getNRows(),BRHS.getNCols());
    result.putTol(BRHS.getTol());
    for(int j=0; j<BRHS.getNCols(); j++)
      for(int i=0; i<BRHS.getNRows(); i++)
	result(i,j)=BRHS(i,j);
    */
    result = BRHS;
  }
  else{ //copy the transpose of B into work matrix "result"
    result.newSize(BRHS.getNCols(),BRHS.getNRows());
    result.putTol(BRHS.getTol());
    for(int i=0; i<BRHS.getNRows(); i++)
      for(int j=0; j<BRHS.getNCols(); j++)
	result(j,i)=BRHS(i,j);    
  }
  int nrhs= static_cast<int>(result.getNCols());
  
  //printf("solve_after_LU_fact: n=%d nrhs=%d\n",n,nrhs);
  int info=0;
  DGETRS_F77(&transA, &n, &nrhs, ALU.ptr(0,0), &lda, ipvt.ptr(0), result.ptr(0,0), &ldb, &info);
  return result;
}


inline void least_squares(MtxDbl& A, MtxDbl& x, MtxDbl& b)
{
  ///\todo Change input b to be pass by value.  It gets overwritten
  /// inside of dgels, and that has caused some problems.  It's probably
  /// better to just copy the whole vector and pass the copy to dgels.
  // Rows in A must == size of b
#ifdef __SURFMAT_ERR_CHECK__
  assert(A.getNRows() == b.getNElems()); 
  // System must be square or over-constrained
  assert(A.getNRows() >= A.getNCols());
#endif
  int n_rows = static_cast<int>(A.getNRows());
  int n_cols = static_cast<int>(A.getNCols());
  // Client may supply a "blank" initialized vector for x
  int lwork = n_rows*n_cols * 2;
  MtxDbl work(lwork);
  // values must be passed by reference to Fortran, so variables must be 
  // declared for info, nrhs, trans
  int info;
  int nrhs=1;
  char trans = 'N';
  x = b; //preserve the original b
  DGELS_F77(&trans,&n_rows,&n_cols,&nrhs,A.ptr(0,0),&n_rows,x.ptr(0),
	    &n_rows,work.ptr(0),&lwork,&info);
  x.reshape(n_cols);
}


inline void least_squares_with_equality_constraints(MtxDbl& A, 
     MtxDbl& x, MtxDbl& c, MtxDbl& B, MtxDbl& d)
{
  int m = static_cast<int>(A.getNRows());
  int n = static_cast<int>(A.getNCols());
  int p = static_cast<int>(B.getNRows());
  int B_cols = static_cast<int>(B.getNCols());
#ifdef __SURFMAT_ERR_CHECK__
  assert(B_cols == n);
  assert(p <= n);
  assert(n <= p + m);
  assert(x.getNElems() == n);
#endif
  int lwork = (m + n + p);
  lwork *= lwork;
  ///\todo Compute optimal blocksize before running dgglse
  MtxDbl work(lwork);
  int info = 0;
  DGGLSE_F77(&m,&n,&p,A.ptr(0,0),&m,B.ptr(0,0),&p,c.ptr(0),d.ptr(0),x.ptr(0),work.ptr(0),&lwork,
	     &info);
#ifdef __SURFMAT_ERR_CHECK__
  if (info != 0) {
    printf("least_squares_with_equality_constraints: Error encountered in DGGLSE\n");
    assert(info == 0);
  }
#endif

  //if (info != 0) throw string("Error in dgglse\n");
}

///finds the eigenvalues and optionally (by default) eigenvectors of a real symmetric matrix, returns a reference to the vector of eigenvalues, this function wraps the LAPACK subroutine DSYEV
inline MtxDbl& eig_sym(MtxDbl& eigvect, MtxDbl& eigval, const MtxDbl& A, char jobz='V') {
  char uplo = 'U';
  int lda=A.getNRows();
#ifdef __SURFMAT_ERR_CHECK__
  assert((0<lda)&&(lda==A.getNCols()));
#endif
  eigvect.copy(A);
  eigval.newSize(lda);eigval.zero();
  int n=lda;
  int info;
  int lwork=-1;
  double work_size;
  DSYEV_F77(&jobz, &uplo, &n, eigvect.ptr(0), &lda, eigval.ptr(0), 
	    &work_size, &lwork, &info);
#ifdef __SURFMAT_ERR_CHECK__
  assert(info==0);
#endif
  lwork=static_cast<int>(work_size);
  MtxDbl work(lwork);
  DSYEV_F77(&jobz, &uplo, &n, eigvect.ptr(0), &lda, eigval.ptr(0), 
	    work.ptr(0), &lwork, &info);

  return eigval;
}

///finds the eigenvalues (but not eigenvectors) of a real symmetric matrix, returns a reference to the vector of eigenvalues, this function wraps the function inline MtxDbl& eig_sym(MtxDbl& eigvect, MtxDbl& eigval, const MtxDbl& A, char jobz='V') which inturn wraps LAPACK subroutine DSYEV
inline MtxDbl& eig_sym(MtxDbl& eigval, const MtxDbl& A) {
  MtxDbl eigvect;
  return (eig_sym(eigvect,eigval,A,'N'));
}

} // end namespace nkm

#endif
