c _______________________________________________________________________
c
c DAKOTA: Design Analysis Kit for Optimization and Terascale Applications
c Copyright (c) 2001, Sandia National Laboratories.
c This software is distributed under the GNU General Public License.
c For more information, see the README file in the top Dakota directory.
c _______________________________________________________________________ 
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c234567
c
c written by Tony Giunta, 26 May 1997,
c
c revised: 11 February 2000
c (1) to interface with DAKOTA KrigingSurf.H/C codes,
c (2) to allow memory management within C++ DAKOTA codes
c        
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c

      subroutine calldot(numdv,numsamp,numnew,
     &  dot_nrwk,dot_nriwk,dot_param_array,dot_single_array,
     &  iflag,theta,thetaLoBnds,thetaUpBnds,xmat,yvect,
     &  xnew,ynew,betahat,RHSterms,MLE,
     &  g,wk,iwk,rprm,iprm,
     &  ipvt,cormat,invmat,Fvect,FRinv,yfb_vect,yfbRinv,r_xhat,work,
     &  work4,iwork,numsamp4)


c    
c variable declarations
c
      implicit double precision (a-h,o-z)

      integer numdv,numsamp,numnew,dot_nrwk,dot_nriwk,dot_param_array,
     &  dot_single_array,iwk(dot_nriwk),ipvt(numsamp),
     &  iprm(dot_param_array), 
     &  nrwk,nriwk,ncon

      double precision wk(dot_nrwk),rprm(dot_param_array),
     &  g(dot_single_array)


c
c declaration of Kriging code variables
c

      integer numsamp4,iwork(numsamp)

      double precision xmat(numsamp,numdv),yvect(1,numsamp),
     &  RHSterms(1,numsamp), betahat, 
     &  MLE, xnew(numnew,numdv), ynew(numnew),
     &  theta(numdv),thetaLoBnds(numdv),thetaUpBnds(numdv)

      double precision cormat(numsamp,numsamp),invmat(numsamp,numsamp),
     &  Fvect(numsamp,dot_single_array),FRinv(numsamp,dot_single_array),
     &  yfb_vect(numsamp,dot_single_array),
     &  yfbRinv(numsamp,dot_single_array),
     &  r_xhat(numsamp,dot_single_array),work(numsamp),
     &  work4(numsamp4)


c
c initialize wk and iwk arrays
c
      nrwk  = dot_nrwk
      nriwk = dot_nriwk
      do 200 i = 1,nrwk
        wk(i) = 0.0
  200 continue
      do 210 i = 1,nriwk
        iwk(i) = 0
  210 continue


c
c initialize rprm and iprm arrays
c
      do 10 i=1,dot_param_array
        rprm(i)=0.0
10      iprm(i)=0


c
c set convergence tolerance criteria to other than default
c
      iprm(4) = 10 
      iprm(9) =  6 
      rprm(1) = -1.0d-1
      rprm(2) =  1.0d-3
      rprm(4) =  1.0d-4


c define method,numdv,ncon
c
c Unconstrained: 
c 0=BFGS
c 1=FR
c
c Constrained:
c 0,1 = modified method of feasible directions
c 2   = sequential linear program
c 3   = sequential quadratic program
C
      ncon   = 0
      method = 0

c print flag
      iprint=4

c set flag: minimize (minmax=0 or -1), maximize (minmax=1)
      minmax=1

c initialize info to zero
      info=0
      
c initial value of x, along with lower and upper bounds

c      do 90 i=1,numdv
c         x(i)  = theta(i)
c         xl(i) = thetaLoBnds(i)
c         xu(i) = thetaUpBnds(i)
c 90   continue

CCC
CCC 'ngotoz' was added to DOT by Vanderplaats, et al, to provide
CCC a better interface with DAKOTA (I'm not sure what info is passed
CCC withing ngotoz)
CCC
      ngotoz = 0

c
c call DOT to start optimization
c
 100  call dot(info,ngotoz,method,iprint,numdv,ncon,
     &   theta,thetaLoBnds,thetaUpBnds,
     &   obj,minmax,g,rprm,iprm,wk,nrwk,iwk,nriwk)

c
c compute objective function for given design variables
c
      iprm(1) = 0
      iprm(2) = 0
      iprm(3) = 0
      iprm(5) = 0
      iprm(10) = 0

      if( info .ne. 0 ) then
c        do 110 i = 1,numdv
c           theta(i) = x(i)
c110     continue
         call krigmodel(numdv,numsamp,numnew,iflag,theta,xmat,
     &     yvect,xnew,ynew,betahat,RHSterms,MLE,
     &     ipvt,cormat,invmat,Fvect,FRinv,yfb_vect,yfbRinv,r_xhat,work,
     &     work4,iwork,numsamp4)
         obj = MLE
         go to 100
      endif

c
c info=0 is a flag that the optimizer is done
c copy the design variables and objective function value back to
c the original variable names
c
c     do 120 i = 1,numdv
c        theta(i) = x(i)
c120  continue

c
c make one final call to krigmodel to compute all of the
c needed kriging terms using the final (optimal) values of the
c correlation parameters
c
      call krigmodel(numdv,numsamp,numnew,iflag,theta,xmat,
     &  yvect,xnew,ynew,betahat,RHSterms,MLE,
     &  ipvt,cormat,invmat,Fvect,FRinv,yfb_vect,yfbRinv,r_xhat,work,
     &  work4,iwork,numsamp4)
      MLE = obj

c
c return to calling program
c         
      return
      end
