c _______________________________________________________________________
c
c DAKOTA: Design Analysis Kit for Optimization and Terascale Applications
c Copyright (c) 2001, Sandia National Laboratories.
c This software is distributed under the GNU General Public License.
c For more information, see the README file in the top Dakota directory.
c _______________________________________________________________________ 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c234567
*
      subroutine krigmodel(numdv,numsamp,numnew,
     &   iflag, theta, xmat, yvect, xnew, ynew, betahat, RHSterms, MLE,
     &   ipvt, cormat, invmat, Fvect, FRinv, yfb_vect, yfbRinv,
     &   r_xhat, work, work4, iwork, numsamp4)

*
* This subroutine acts as an interface between a numerical optimizer by
* providing the MLE estimate of the DACE correlation parameter 'theta'
* based on the data contained in the input file opened below.
*
* This routine calls the correlation subroutine (krigmodelbuild) to 
* evaluate the correlation matrix and kriging model parameters (iflag=1).
* Also, this routine evaluates the kriging model at a list of new sample
* sites and returns f_hat(x_new_samp) (iflag=2)
*
* Author: Tony Giunta, 12 May 1997
*                      13 Sept 1999, rev - make 'theta' a vector, remove 'ibeta'
*
******************************************************************************
*
* Input Variables:
* ----------------
*   xnew    = vector of length 'numdv' where the DACE model is to be evaluated
*
* Output Variables:
* -----------------
*   MLE    = the Maximum Likelihood Estimate (MLE) of 'theta,' -- the 
*               correlation parameter.
*
* Parameter Variables (to be set by user at run time):
* ----------------------------------------------------
*   numdv   = number of variables
*   numsamp = number of data samples from which the correlation matrix and
*             the DACE model parameters are calculated
*   numnew  = number of sample site where the unknown response value will be
*             calculated (always = 1 in hsct code)
*   theta   = user-supplied value of the correlation parameter 'theta' 
*
* Local Variables:
* ----------------
*   DOUBLE PRECISION
*   ----------------
*   xmat     = numdv x numsamp of sample site locations 
*   cormat   = correlation matrix (numsamp x numsamp)
*   invmat   = inverse of the correlation matrix (numsamp x numsamp)
*   Fvect    = matrix (1 x numsamp) of constant terms (all = 1 in 'correlate')
*   FRinv    = matrix product of 'Fvect' and 'invmat'
*   yvect    = matrix (1 x numsamp) of response values
*   yfb_vect = matrix (1 x numsamp) resulting from ('yvect' - 'Fvect'*'betahat')
*   yfbRinv  = matrix (1 x numsamp) resulting from 
*                 ('yvect' - 'Fvect'*'betahat')*'invmat'
*   RHSterms = matrix product of 'invmat' and 'yfb_vect'
*   r_xhat   = matrix (1 x numsamp) created by using the vector 'xnew' in the
*                 correlation function
*   betahat  = estimate of the constant term in the DACE model (beta)
*   sigmahat = estimate of the variance (sigma) term in the data
*   MLE      = scalar valued function determined by choice of 'theta' 
*                 used for comparing to maximum MLE value found from prior 
*                 analysis of the response data
*   work     = vector of length 'numsamp' used as temporary storage by the
*                 LAPACK subroutine DGEDI
*
*   INTEGER
*   -------
*   iflag    = flag variable:
*                 = 1 --> calculate kriging modeling terms & correlation matrix
*                 = 2 --> calculate predicted response with the kirging model
*			  using previously calculated DACE model parameters 
*                         found when iflag=1
*   ipvt     = vector of length 'numsamp' of pivot locations used in LAPACK
*                 subroutines DGEDI and DGEFA
*
*
* Notes:
*   1. User must supply variable values in parameter statement prior to 
*      executing this code.
*   2. User must supply the data from which the correlation matrix is
*      calculated in the matrix form (numdv + 1 columns by numsamp rows).
*    
*   x(1,1)      x(1,2)          ...  x(1,numdv)        response_value(1)
*      |            |           ...     |                    |
*      |            |           ...     |                    |
*   x(numsamp,1)  x(numsamp,2)  ...  x(numsamp,numdv)  response_value(numsamp)
*
*
****************************************************************************
*

      integer i,j,numdv,numsamp,numsamp4,ipvt(numsamp),iflag,
     &   iwork(numsamp)

      double precision xmat(numsamp,numdv),cormat(numsamp,numsamp),
     &   invmat(numsamp,numsamp),Fvect(1,numsamp),FRinv(1,numsamp),
     &   yvect(1,numsamp),yfb_vect(1,numsamp),yfbRinv(1,numsamp),
     &   RHSterms(1,numsamp),r_xhat(1,numsamp), betahat, sigmahat, 
     &   MLE, work(numsamp), xnew(numnew,numdv), ynew(numnew),
     &   theta(numdv),work4(numsamp4)

     
*
* initialize the DACE modeling parameters if needed (iflag=1)
*
      if( iflag .eq. 1 ) then
         betahat  = 0.0d0
         MLE      = 0.0d0
         sigmahat = 0.0d0
      endif

*
* initialize y_new(x)
*
      do 130 i = 1,numnew
         ynew(i) = 0.0d0
  130 continue

*
* initialize all other vectors and matrices
*
      do 200 i=1,numsamp
         do 210 j = 1,numsamp
            cormat(i,j) = 0.0d0
            invmat(i,j) = 0.0d0
  210    continue      
  200 continue 

      do 220 i=1,numsamp
         Fvect(1,i)    = 0.0d0
         FRinv(1,i)    = 0.0d0
         yfb_vect(1,i) = 0.0d0
         yfbRinv(1,i)  = 0.0d0 
         r_xhat(1,i)   = 0.0d0
         work(i)       = 0.0d0         
         if( iflag .eq. 1 ) RHSterms(1,i) = 0.0d0
         ipvt(i)       = 0
  220 continue          

*
* call subroutine to calculate the inverse of the correlation matrix 
* and the correlation parameters
*
* iflag = 1, for the given 'theta' find the kriging model parameters
*
* iflag = 2, evaluate Y(x_new), given x_new, correlation parameters,
* and the existing sample site data (xmat)
*
      if( iflag .eq. 1 ) then
      call krigmodelbuild (xmat,cormat,
     &   invmat,Fvect,FRinv,yvect,yfb_vect,yfbRinv,RHSterms,r_xhat,
     &   numsamp, numdv, numnew,betahat,sigmahat,MLE, work,theta,ipvt,
     &   xnew,ynew,work4,iwork,numsamp4)   
      endif


      if( iflag .eq. 2 ) then
      call krigmodeleval(xnew,xmat,r_xhat,betahat,RHSterms,numsamp,
     &   numdv,numnew,theta,ynew)      
      endif

      return
      end


***********************************************************************
c234567

      subroutine krigmodelbuild (xmat, cormat,
     &   invmat, Fvect, FRinv, yvect, yfb_vect, yfbRinv,
     &   RHSterms, r_xhat, numsamp, numdv, numnew,
     &   betahat, sigmahat, MLE, work, theta, ipvt, xnew, ynew,
     &   work4, iwork, numsamp4)
*     
*     
* This subroutine calculates the kriging correlation matrix 
* using an exponential correlation function (see Boeing/IBM/Rice Report)
* along with a user-supplied value of the correlation parameter 'theta.'
*
* Tony Giunta, 12 May 1997
*              13 Sept. 1999, rev - make 'theta' a vector, remove 'ibeta'
* 
*
***************************************************************************
*
* Needed External Files: 
* ----------------------
*   LAPACK subroutines DGETRF, DGETRI, DGECON (FORTRAN77)
*   BLAS subriutine DGEMM (FORTRAN77)
*
* Inputs: (variables defined in calling subroutine )
* -------
*   xmat    
*   cormat
*   Fvect
*   FRinv
*   yvect  
*   yfb_vect
*   yfbRinv
*   RHSterms
*   r_xhat
*   numsamp 
*   numdv   
*   numnew
*   sigmahat
*   MLE
*   work
*   ipvt
*   xnew
*   ynew
*
*
* Outputs: (variables defined in calling subroutine )
* --------
*   invmat  
*   betahat 
*   theta   
*
*
******************************************************************************
*
      character*1 norm

      integer numdv,numsamp,numnew,i,j,k,ipvt(numsamp),info

      integer numsamp4,iwork(numsamp)

      double precision xmat(numsamp,numdv),
     &   cormat(numsamp,numsamp),
     &   invmat(numsamp,numsamp),Fvect(1,numsamp),FRinv(1,numsamp),
     &   yvect(1,numsamp),beta_num,beta_den,betahat,sigmahat,MLE,
     &   yfb_vect(1,numsamp),yfbRinv(1,numsamp),
     &   work(numsamp),RHSterms(1,numsamp),
     &   xnew(numnew,numdv), r_xhat(1,numsamp), theta(numdv),
     &   detR,sum,ynew(numnew)

      double precision anorm,rcond,cond_num
      double precision work4(numsamp4)

*
* calculate the correlation matrix, its inverse, and correlation parameters
*      
*
* calculate terms in the correlation matrix
*      
 
       do 200 i = 1,numsamp
          do 210 j = i,numsamp
             if( i .eq. j ) then
                cormat(i,j) = 1.0d0
                invmat(i,j) = 1.0d0
             else
                sum = 0.0d0
                do 220 k = 1,numdv
                   sum = sum + theta(k)*
     &                   (xmat(i,k)-xmat(j,k))*(xmat(i,k)-xmat(j,k))
  220           continue     
                cormat(i,j) = exp( -1.0d0*sum )
                cormat(j,i) = cormat(i,j)
                invmat(i,j) = cormat(i,j)
                invmat(j,i) = cormat(i,j)
             endif
  210     continue     
  200  continue


*
* for the correlation matrix calculate the following:
*    1. condition number
*    2. determinant
*    3. inverse
* 
* note: all of these require an LU factorization of 'cormat'
*

*
* perform the LU factorization of the correlation matrix
*
         call dgetrf( numsamp, numsamp, invmat, numsamp, ipvt, info )

*
* compute the condition number, but must first compute the 1-norm
* of the correlation matrix
*
* note: 1-norm = sum_over_i sum_over_j abs(a_ij)
*
         anorm = 0.0d0
         do 222 i = 1,numsamp
            do 224 j = 1,numsamp
               anorm = anorm + dabs( cormat(i,j) )
 224        continue     
 222     continue



* set norm flag:  '0' or '1' ==> 1-norm, 'I' ==> infinity norm
         norm = '1'
         call dgecon( norm, numsamp, invmat, numsamp, anorm, rcond, 
     &                work4, iwork, info )
         cond_num = 1.0/rcond  
      if( cond_num .ge. 1.0d+32 ) then
         write(*,*)
         write(*,*)"***************************************************"
         write(*,*)"Error in Kriging Model: Ill-conditioned Corr.Matrix"
         write(*,*)"***************************************************"
         write(*,*)
         stop
      endif
  
*
* calculate the determinant of the correlation matrix using the
* LU factorization from LAPACK routine DGETRF:
*    det(cormat) = det(LU) = det(L)*det(U)
*                = 1*Product_over_i(U_ii)   (i.e., diag elements of U)
*
       detR = 0.0d0
       do 230 i = 1,numsamp
          if( i .eq. 1 ) then
             detR = invmat(i,i)
          else
             detR = detR*invmat(i,i)
          endif
 230   continue
       detR = dabs(detR)

c      if( detR .lt. 1.0d-30 ) then
c         write(*,*)
c         write(*,*)"**************************************************"
c         write(*,*)"Error in Kriging Model: Determinant[Corr. Matrix]"
c         write(*,*)"Det[R] = ",detR
c         write(*,*)"**************************************************"
c         write(*,*)
c         stop
c      endif

*
* compute the inverse of the correlation matrix using LAPACK 
* routine DGETRI
*
      info = 0
      call dgetri (numsamp,invmat,numsamp,ipvt,work,numsamp,info)

*
* calculate the terms needed to determine the MLE value
*
      do 300 i = 1,numsamp
         Fvect(1,i) = 1.0d0
  300 continue          
*
      call dgemm('n','n',1,numsamp,numsamp,1.0d0,Fvect,1,invmat,
     &   numsamp,0.0d0,FRinv,1)
      call dgemm('n','n',1,1,numsamp,1.0d0,FRinv,1,Fvect,numsamp,
     &   0.0d0,beta_den,1)
      call dgemm('n','n',1,1,numsamp,1.0d0,FRinv,1,yvect,numsamp,
     &   0.0d0,beta_num,1)

      if( abs(beta_den) .ge. 1.0d-6 ) then
         betahat = beta_num / beta_den
      else
         write(*,*)
         write(*,*)"**************************************************"
         write(*,*)"    Error in Kriging Model: Beta_hat Estimate"
         write(*,*)"**************************************************"
         write(*,*)
         stop
      endif

cccccccccccccccccccccccccccc
c
c Debug code
c
c      if( abs(beta_den) .ge. 1.0d-6 ) then
c         betahat = beta_num / beta_den
c      else
c         write(*,*)
c         write(*,*) 'Error in calculating Beta_hat in kriging mode.'
c         write(*,*)
c         write(*,*) 'beta_num and beta_den are', beta_num,beta_den
c         write(*,*) 'condition number = ', cond_num
c         write(*,*) 'det(R) = ', detR
c         do 301 i=1,numdv
c            write(*,*)'theta(i) = ', theta(i)
c 301     continue
c         write(*,*)
c         do 302 i=1,numsamp
c            write(*,*)'i, x(i,1_), x(i,2) = ', i,xmat(i,1),xmat(i,2)
c 302     continue
c         write(*,*)
c         do 303 i=1,numsamp
c            write(*,*) (cormat(i,j),j=i,numsamp)
c 303     continue
c         stop
c      endif
c      
ccccccccccccccccccccccc


      do 310 i = 1,numsamp
         yfb_vect(1,i) = yvect(1,i) - betahat*Fvect(1,i)
  310 continue
      
      call dgemm('n','n',1,numsamp,numsamp,1.0d0,yfb_vect,1,invmat,
     &   numsamp,0.0d0,yfbRinv,1)
      call dgemm('n','n',1,1,numsamp,1.0d0,yfbRinv,1,yfb_vect,
     &   numsamp,0.0d0,sigmahat,1)
      sigmahat = sigmahat / numsamp

      sigmahat = dabs(sigmahat)

      if( sigmahat .lt. 1.0d-6 ) then
         write(*,*)
         write(*,*)"************************************************"
         write(*,*)"  Warning: Kriging model sigma_hat term is zero "
         write(*,*)"     and the MLE value is also zero.            "
         write(*,*)"  This can happen if all data points have the   "
         write(*,*)"     same response value.                       "
         write(*,*)"************************************************"
         write(*,*)
         MLE = 0.0d0
      else
         MLE = -0.5d0*(numsamp*log(sigmahat) + log(detR))
      endif





* 
* calculate the "right hand side" terms in the kriging model
* i.e., R_inverse * (y-f*beta),  and return this quantity to
* calling program
*
      call dgemm('n','n',numsamp,1,numsamp,1.0d0,invmat,numsamp,
     &   yfb_vect,numsamp,0.0d0,RHSterms,numsamp)
      

      
      return     
      end      


**************************************************************************
*
      subroutine krigmodeleval(xnew,xmat,r_xhat,betahat,RHSterms,
     &   numsamp,numdv,numnew,theta,ynew)
*     
*     
* Use kriging interpolating model to predict response values at unsampled
* locations
*
* Tony Giunta, 12 May 1997
*              13 Sept 1997, rev - make 'theta' a vector, remove 'ibeta'
*
**************************************************************************
*
* Inputs: (variables defined in calling subroutine)
* -------
*   xnew
*   xmat
*   r_xhat
*   betahat
*   RHSterms
*   numsamp
*   numdv
*   theta
*
* Outputs: (variables defined in calling subroutine)
* --------
*   ynew
*
* Local Variables:
* ----------------
*   sum   = temporary variable used for calculating the terms in the
*              vector 'r_xhat'
*   yeval = scalar value resulting from matrix multiplication of
*              'r_xhat' * 'RHSterms'
*
*************************************************************************
*
      double precision xnew(numnew,numdv),r_xhat(1,numsamp),
     &   xmat(numsamp,numdv),RHSterms(1,numsamp),betahat,
     &   theta(numdv),sum,ynew(numnew),yeval
      
      integer i,j,k,numdv,numsamp,numnew
*
* calculate the vector r(x)
*     
      do 200 i = 1,numnew
      
         do 110 j = 1,numsamp
            sum = 0.0d0
            do 120 k = 1,numdv
               sum = sum + theta(k)*
     &            (xnew(i,k)-xmat(j,k))*(xnew(i,k)-xmat(j,k))
  120          continue   
            r_xhat(1,j) = exp( -1.0d0*sum )
  110    continue 
*
* calculate the estimate of Y, i.e., Y_hat(x) using the previously
* found kriging model terms 'betahat' and 'RHSterms'.
*
         yeval = 0.0d0
         call dgemm('n','n',1,1,numsamp,1.0d0,r_xhat,1,RHSterms,
     &      numsamp,0.0d0,yeval,1)
      
         ynew(i) = yeval + betahat

  200 continue

      return
      end


      

