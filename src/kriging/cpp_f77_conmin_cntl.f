c _______________________________________________________________________
c
c DAKOTA: Design Analysis Kit for Optimization and Terascale Applications
c Copyright (c) 2001, Sandia National Laboratories.
c This software is distributed under the GNU General Public License.
c For more information, see the README file in the top Dakota directory.
c _______________________________________________________________________ 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c234567
c
c Author: Tony Giunta, October 2001
c
c
      subroutine callconmin(
     &  thetaVars,
     &  thetaLoBnds, 
     &  thetaUpBnds,
     &  g,
     &  SCAL,DF,A,S,G1,G2,B,C,ISC,IC,MS1,N1,N2,N3,N4,N5,
     &  DELFUN,DABFUN,FDCH,FDCHM,CT,CTMIN,
     &  CTL,CTLMIN,ALPHAX,ABOBJ1,THETA,
     &  MLE,
     &  numdv,numcon,NSIDE,IPRINT,NFDG,NSCAL,LINOBJ,ITMAX,
     &  ITRM,ICNDIR,IGOTO,NAC,conminInfo,INFOG,ITER,
     &  numsamp,numnew,iflag,xmat,
     &  yvect,xnew,ynew,betahat,RHSterms,
     &  ipvt,cormat,invmat,Fvect,FRinv,yfb_vect,yfbRinv,r_xhat,work,
     &  work4,iwork,numsamp4,one_element)

c    
c CONMIN variables and shared variables
c
      implicit double precision (a-h,o-z)

      integer numdv,numsamp,numnew,numcon,one_element,N1,N2,N3,N4,N5,
     &  NSIDE,IPRINT,NFDG,NSCAL,LINOBJ,ITMAX,ITRM,ICNDIR,IGOTO,NAC,
     &  INFOG,ITER,conminInfo,iflag

      integer MS1(N5),ISC(N2),IC(N3),ipvt(numsamp)

      double precision obj,DELFUN,DABFUN,FDCH,FDCHM,CT,CTMIN,
     &  CTL,CTLMIN,ALPHAX,ABOBJ1,THETA

      double precision S(N1),G1(N2),G2(N2),B(N3*N3),C(N4),SCAL(N1),
     &  DF(N1),A(N1*N3)

c
c kriging code variables
c

      integer numsamp4,iwork(numsamp)

      double precision xmat(numsamp,numdv),yvect(1,numsamp),
     &  RHSterms(1,numsamp), betahat, 
     &  MLE, xnew(numnew,numdv), ynew(numnew),
     &  thetaVars(N1),thetaLoBnds(N1),thetaUpBnds(N1),
     &  g(N2)

      double precision cormat(numsamp,numsamp),invmat(numsamp,numsamp),
     &  Fvect(numsamp,one_element),
     &  FRinv(numsamp,one_element),
     &  yfb_vect(numsamp,one_element),
     &  yfbRinv(numsamp,one_element),
     &  r_xhat(numsamp,one_element),work(numsamp),
     &  work4(numsamp4)

c
c call CONMIN to start optimization
c
CCC 100  call dot(info,ngotoz,method,iprint,numdv,ncon,
CCC     &   theta,thetaLoBnds,thetaUpBnds,
CCC     &   obj,minmax,g,rprm,iprm,wk,nrwk,iwk,nriwk)

      OBJ = MLE
      

 100  call conmin(
     &      thetaVars,thetaLoBnds,thetaUpBnds,
     &      g,
     &      SCAL,DF,A,S,G1,G2,B,C,
     &      ISC,IC,MS1,N1,N2,N3,N4,N5,
     &      DELFUN,DABFUN,FDCH,FDCHM,CT,CTMIN,
     &      CTL,CTLMIN,ALPHAX,ABOBJ1,THETA,
     &      OBJ,
     &      numdv,numcon,NSIDE,IPRINT,NFDG,NSCAL,LINOBJ,ITMAX,
     &      ITRM,ICNDIR,IGOTO,NAC,conminInfo,INFOG,ITER)

c
c IGOTO = 0 is a flag that CONMIN has finished
c
      if (IGOTO .ne. 0) then
c
c compute objective function value for current design variables
c
         if( conminInfo .ne. 0 ) then
            call krigmodel(numdv,numsamp,numnew,iflag,thetaVars,xmat,
     &           yvect,xnew,ynew,betahat,RHSterms,MLE,
     &           ipvt,cormat,invmat,Fvect,FRinv,yfb_vect,yfbRinv,
     &           r_xhat,work,work4,iwork,numsamp4)
            OBJ = -1.0d0*MLE
            go to 100
         endif
      endif

c
c make one final call to krigmodel to compute all of the
c needed kriging terms using the final (optimal) values of the
c correlation parameters
c
      call krigmodel(numdv,numsamp,numnew,iflag,thetaVars,xmat,
     &  yvect,xnew,ynew,betahat,RHSterms,MLE,
     &  ipvt,cormat,invmat,Fvect,FRinv,yfb_vect,yfbRinv,r_xhat,work,
     &  work4,iwork,numsamp4)
      MLE = -1.0d0*OBJ

c
c return to calling program
c         
      return
      end
