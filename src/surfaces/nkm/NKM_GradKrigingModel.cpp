#include "NKM_SurfPack.hpp"
#include "NKM_GradKrigingModel.hpp"
//#include "Accel.hpp"
#include "NKM_LinearRegressionModel.hpp"
#include <cmath>
#include <iostream>
#include <cfloat>
//#define __TIME_PIVOT_CHOLESKY__
#ifdef __TIME_PIVOT_CHOLESKY__
#include <sys/time.h>
#endif
namespace nkm {

using std::cout;
using std::cerr;
using std::endl;
using std::ostringstream;

//#define __DEBUG_MY_CHOL__
//#define __DEBUG_MY_CHOL2__
//#define __NKM_UNBIASED_LIKE__
#define __PIVOT_CHOL_KRIGING_R__
//#define __PIVOT_CHOL_GEK_R__

//#define __GRADKRIGING_DER_TEST__


void GradKrigingModel::getBaseIEqnKeep() {
  //get the index to the equations that we are required to keep, i.e. the equations that we are required to match exactly.  These are the function values at all points, PLUS the first derivatives at "Anchor Points"

  iOrderEqnTest.newSize(numEqnAvail,1);
  iPointOrderTest.newSize(numPoints,1);
  int npts_so_far=0;
  if(ifHaveAnchorPoint==true) {
    iPointOrderTest(0,0)=iAnchorPoint;
    npts_so_far=1;
  }
  else
    iAnchorPoint=-1;
  //if(ifHaveAnchorPoint==false)
  //std::cerr << "ifHaveAnchorPoint==false" <<std::endl;

  //std::cerr << "C: ifHaveAnchorPoint=" << ifHaveAnchorPoint <<
  //" iAnchorPoint=" << iAnchorPoint << std::endl;
  
  //printf("ifHaveAnchorPoint=%d iAnchorPoint=%d\n",
  // ifHaveAnchorPoint, iAnchorPoint);
  fflush(stderr);

  for(int i=0; i<numPoints; ++i) 
    if(i!=iAnchorPoint) {
      iPointOrderTest(npts_so_far,0)=i;
      ++npts_so_far;
    }
  
  iEqnKeep.newSize(numEqnAvail,1);
  iptIderKeep.newSize(numEqnAvail,2);
}

#if defined(__PIVOT_CHOL_KRIGING_R__)
/// select which derivative equations to add based on whether that information is already present in R, if the information is already present in R then adding the equation will make R poorly conditioned which means we can use rcond(R) as a criteria whether to accept or reject equations.  This function uses pivoted Cholesky to reorder points in Kriging's (not GEK's) R matrix so that the points with the most "new" information come first and the points that contain the least new information come last.  This reordering is then applied to WHOLE POINTS (a whole point is a function value at a point immediately followed by the derivatives at the same point) in GEK's R, LAPACK's highly efficient LEVEL3 BLAS non-pivoted Cholesky is used on the reordered GEK R matrix, and equations are dropped of the end of R until it is not "too ill-conditioned" (until 2^-40<rcond(GEK R)). A bisection search, at most ceil(log2(# of equations)) evaluations of rcond, is then used to efficiently find the last equation in R that can be retained (this is all based of the same Cholesky factorization of the reordered GEK R but uses a different number of equations to compute rcond).  Since the points/equations dropped are the ones with the least new information they should be the ones that when dropped are both the easiest to predict AND the ones that most improve the conditioning of R. The pivoted Cholesky implementation used is located in NKM_PivotChol.f, it's a modification of the implementation found in C. Lucas, "LAPACK-style Codes for Level2 and 3 Pivoted Cholesky Factorizations", Numerical Analysis Report No. 442, February 2004, from the Manchester Center for Computational Mathematics. I KRD downloaded it from http://www.maths.manchester.ac.uk/~nareports/narep442.pdf. C Lucas's implementation was the best pivoting Cholesky algorithm known to David Day.  I modified it to improve efficiency, in tests the modified version was between 5% and 10% faster than Lucas's level 2 (and level 3) implementation and I (KRD) believe that Lucas's Level 3 implementation was ALWAYS defaulting to his level 2 implementation.  Note that perfoming pivoted Cholesky on Kriging's R is about (1+numVarsr)^3 (assuming that a function value plus gradient is used for GEK) than doing so for GEK's R __AND__ it produces more accurate GEK emulators.  Doing pivoted Cholesky on GEK's R tends to sort function values to the end, which get dropped, while some derivative equations at the same points are retained (that had a detrimental effect on emulator prediction accuracy)
void GradKrigingModel::equationSelectingPrecondCholR(){ 
#ifdef __TIME_PIVOT_CHOLESKY__
  ++n_pivot_cholesky_calls;
  struct timeval tv;
  gettimeofday(&tv, NULL); 
  long int start_sec_block4, start_usec_block4;
  long int start_sec_cholesky=tv.tv_sec;
  long int start_usec_cholesky=tv.tv_usec;
  long int stop_sec_cholesky, stop_usec_cholesky;
  long int start_sec_rcond, start_usec_rcond;
  long int stop_sec_rcond, stop_usec_rcond;
#endif  
  double min_allowed_rcond=1.0/maxCondNum;
  //double min_allowed_pivot_est_rcond=256.0/maxCondNum;

  scaleRChol.newSize(numEqnAvail,2); 
  for(int i=0; i<numEqnAvail; ++i) {
    scaleRChol(i,1)=std::sqrt(R(i,i));
    scaleRChol(i,0)=1.0/scaleRChol(i,1);
  }
  //precondition R
  for(int j=0; j<numEqnAvail; ++j) {
    for(int i=0; i<numEqnAvail; ++i)
      R(i,j)*=scaleRChol(i,0)*scaleRChol(j,0);
    R(j,j)=1.0; //round off error on the diagonal can screw up pivoting 
    //cholesky... since there is zero correlation between an individual 
    //point's function value and its derivatives and since there are 
    //only ones on the diagonal of the preconditioned R, then pivoting
    //cholesky will not pivot the first point's function value or 
    //derivatives (because there won't be any diagonal entries larger 
    //than one, and the preconditioned correlation between first point's 
    //function value and derivatives is the identity matrix) round off
    //error on the diagonal can screw with that and I rely on those 
    //assumptions to ensure that the anchor point's derivatives won't be
    //"dropped"
  }

  RChol.newSize(numEqnAvail,numEqnAvail);
  int ld_RChol=RChol.getNRowsAct();

  //if the user specifies an anchor point we will do an initial 
  //pivoting cholesky with the anchor point's function value and 
  //derivatives first, followed by the function values only at 
  //the rest of the data points
  //
  //if the user doesn't specify an anchor point we will do an 
  //initial pivoting with only the function values
  //
  //the purpose is to pivot/order the points in a good way, i.e.
  //by pivoting Cholesky, based on a small number of equations, 
  //rather than all equations, because pivoting Cholesky is 
  //significantly more expensive than the highly optimized 
  //non-pivoting LAPACK Cholesky which (if needed) will be 
  //applied to the full set of equations after the points have 
  //been ordered according to the previous pivoting.
  //
  //because the function value and its derivatives at an 
  //individual point are uncorrelated with each other and we
  //have normalized R so it has all ones on it's diagonals,
  //pivoting Cholesky WILL NOT PIVOT the first point's 
  //function values or derivatives (the first equation that 
  //can be pivoting is the second point's function value, read
  //as the firs equation that does not have zero correlation
  //with the first point or any of its derivatives)

  int num_eqn_in_grp=numVarsr+1; 
  int num_eqn_test=0;
  lapackRcondR.newSize(numEqnAvail,1);
  lapackRcondR.zero(); 

  //the first point
  int ipt=iPointOrderTest(0,0);
  iOrderEqnTest.newSize(numEqnAvail,1);
  iOrderEqnTest(num_eqn_test,0)=ipt;
  lapackRcondR(num_eqn_test,0)=1.0;
  ++num_eqn_test;

  //if the first point is an anchor point add its derivatives next
  //if(ifHaveAnchorPoint==true) 
  //for(; num_eqn_test<num_eqn_in_grp; ++num_eqn_test) {
  //  iOrderEqnTest(num_eqn_test,0)=num_eqn_test*numPoints+ipt;
  //  lapackRcondR(num_eqn_test,0)=1.0; //identity matrix so far
  //}
  int iprev_lapack_rcondR=num_eqn_test-1; //will use for a bisection
  //search for the last equation we can keep without it being too
  //poorly conditioned 
  
  //then add the rest of the points' function values but not 
  //their derivatives
  for(int i=1; i<numPoints; ++i, ++num_eqn_test) {
    iOrderEqnTest(num_eqn_test,0)=iPointOrderTest(i,0);
    //    if(iOrderEqnTest(num_eqn_test,0)>=numEqnAvail) {
    //printf("numPoints=%d i=%d num_eqn_test=%d iPointOrderTest(i)=%d\n",
    //     numPoints,i,num_eqn_test,iPointOrderTest(i,0));
    //assert(false);
    //}
  }

  //copy the selected equations from the (preconditioned) R
  //matrix to RChol
  for(int j=0; j<num_eqn_test; ++j) {
    int jeqn=iOrderEqnTest(j,0);
    for(int i=0; i<num_eqn_test; ++i) {
      int ieqn=iOrderEqnTest(i,0);
      RChol(i,j)=R(ieqn,jeqn);
    }
  }

  int info=0;
  char uplo='B'; //'B' means we have both halves of R in RChol so the 
  //fortran doesn't have to copy one half to the other (NKM_PivotChol.f was faster than lucas's implentation because having both halves of R allowed more sequential access of memory, which improves cache access performance)
  iPivot.newSize(num_eqn_test,1);  
  numEqnKeep=numEqnAvail;
  PIVOTCHOL_F77(&uplo, &num_eqn_test, RChol.ptr(0,0), &ld_RChol,
    		iPivot.ptr(0,0), &numEqnKeep, &min_allowed_rcond, 
		//&min_allowed_pivot_est_rcond, 
		&info); 
  //printf("Pivoting Cholesky: says num_eqn_test=%d numEqnKeep=%d\n",num_eqn_test,numEqnKeep);

  //for(int i=0; i<numEqnKeep; ++i)
  for(int i=0; i<num_eqn_test; ++i) {
    iPivot(i,0)-=1;
    //printf("iPivot(%d)=%d\n",i,iPivot(i,0));
  }

  oneNormPrecondR.newSize(numEqnAvail,1);
  sumAbsColPrecondR.newSize(numEqnAvail,1); //used in computing the one norm 
  iEqnKeep.newSize(numEqnAvail,1);
  if(false) {
    //if(false) => do not give special treatment to the case where the 
    //__Kriging__ (as opposed to GEK)  R matrix is singular.  Inside this
    //if statement meant that derivative equations would not be used when
    //the __Kriging__ R was singular, test showed that prediction accuracy
    //is better with the default behavior.  The default behavior is to 
    //pivot Kriging R, use the same reordering for WHOLE points (a function 
    //value at a point is immediately followed by the derivatives at the 
    //same point) in GEK R, and then equations at the end of the reordered 
    //GEK R are dropped off until it is not "too ill-conditioned".  
    //The default behavior will not happen if the code enters this if 
    //statement

    //if(numEqnKeep<num_eqn_test) {
    //printf("case 1\n");

    //for this set of correlation lengths, it was so poorly conditioned 
    //it couldn't use all of the function values, let alone all the 
    //derivatives so we are going to, as necessary, drop points of the
    //end until it isn't poorly conditioned

    for(int i=0; i<numEqnKeep; ++i)
      iEqnKeep(i,0)=iOrderEqnTest(iPivot(i,0),0);

    //precompute and store the one norm after adding every equation 
    //in the pivoted Cholesky order, this is prepwork for the 
    //bisection search
    int jeqn=iEqnKeep(0,0);
    for(int i=0; i<numEqnAvail; ++i) {
      sumAbsColPrecondR(i,0)=std::fabs(R(iEqnKeep(i,0),jeqn));
    }
    oneNormPrecondR(0,0)=sumAbsColPrecondR(0,0);
    
    double tempdouble;
    for(int j=1; j<numEqnKeep; ++j) {
      jeqn=iEqnKeep(j,0);
      for(int i=0; i<numEqnAvail; ++i) 	
	sumAbsColPrecondR(i,0)+=std::fabs(R(iEqnKeep(i,0),jeqn));
      tempdouble=sumAbsColPrecondR(0,0);
      for(int i=1; i<=j; ++i) 
	if(tempdouble<sumAbsColPrecondR(i,0))
	  tempdouble=sumAbsColPrecondR(i,0);
      oneNormPrecondR(j,0)=tempdouble;
    }
  }
  else{
    //printf("case 2\n");


    //use the point order from the pivoting cholesky but for the whole 
    //points (function values + derivatives), equations will be ordered a 
    //whole point at a time
    numEqnKeep=numEqnAvail;

    int ipt=iEqnKeep(0,0)=iOrderEqnTest(0,0);
    int ieqn=1;
    int jeqn=1;
    for(; jeqn<num_eqn_in_grp; ++jeqn) {
      iEqnKeep(jeqn,0)=numPoints*jeqn+ipt;
      lapackRcondR(jeqn,0)=1.0; //the preconditioned correlation matrix 
      //between the first point's function value and its derivatives is 
      //the identity matrix, which has condition number = 1 (exactly)
    }
    iprev_lapack_rcondR=jeqn-1;

    //if(ifHaveAnchorPoint==true) 
    //ieqn=jeqn; //inside the commented out (for debug) if

    for(; ieqn<num_eqn_test; ++ieqn) {
      ipt=iOrderEqnTest(iPivot(ieqn,0),0);
      for(int j=0; j<num_eqn_in_grp; ++j, ++jeqn) {
	iEqnKeep(jeqn,0)=numPoints*j+ipt;
	//if(iEqnKeep(jeqn,0)>=numEqnAvail){
	//printf("num_eqn_test=%d ieqn=%d ipt=%d j=%d jeqn=%d iEqnKeep(jeqn)=%d\n",
	// num_eqn_test,ieqn,ipt,j,jeqn,iEqnKeep(jeqn,0));
	//assert(false);
	//}
      }
    }
    //printf("numEqnAvail=%d jeqn=%d\n",numEqnAvail,jeqn);

    //copy from default order of R into test order in RChol
    for(int j=0; j<numEqnAvail; ++j) {
      jeqn=iEqnKeep(j,0);
      for(int i=0; i<numEqnAvail; ++i) {
	ieqn=iEqnKeep(i,0);
	RChol(i,j)=R(ieqn,jeqn);
      }
    }

    //precompute and store the one norm after adding every equation 
    //in the pivoted Cholesky whole-point at a time order (where a 
    //whole-point is the function value and its derivatives at an
    //individual point) this is prepwork for the bisection search
    for(int i=0; i<numEqnAvail; ++i)
      sumAbsColPrecondR(i,0)=std::fabs(RChol(i,0));
    oneNormPrecondR(0,0)=sumAbsColPrecondR(0,0);
    
    double tempdouble;
    for(int j=1; j<numEqnAvail; ++j) {
      for(int i=0; i<numEqnAvail; ++i)
	sumAbsColPrecondR(i,0)+=std::fabs(RChol(i,j));
      tempdouble=sumAbsColPrecondR(0,0);
      for(int i=1; i<=j; ++i) 
	if(tempdouble<sumAbsColPrecondR(i,0))
	  tempdouble=sumAbsColPrecondR(i,0);
      oneNormPrecondR(j,0)=tempdouble;
    }

    //do the (highly optimized) LAPACK Cholesky Decomposition of all 
    //the equations (but sorted into the point order determined by 
    //the pivoting cholesky above) 
    uplo='L';
    DPOTRF_F77(&uplo,&numEqnAvail,RChol.ptr(0,0),&ld_RChol,&info);
  }
  
  //at this point numEqnKeep is not YET the final number of equations 
  //we will keep, only the maximum number of candidates to keep
  //
  //we will use the LAPACK rcond function with a bisection search to
  //find the last equation that isn't too poorly conditioned, but 
  //first we need the rcond at the point the Cholesky decomposition 
  //terminated.

#ifdef __TIME_PIVOT_CHOLESKY__
  gettimeofday(&tv, NULL); 
  start_sec_rcond=tv.tv_sec;
  start_usec_rcond=tv.tv_usec;
#endif

  uplo='L';
  rcondDblWork.newSize(3*ld_RChol,1);
  rcondIntWork.newSize(ld_RChol,1);
  int icurr_lapack_rcondR=numEqnKeep-1;
  DPOCON_F77(&uplo,&numEqnKeep,RChol.ptr(0,0),&ld_RChol,
	     oneNormPrecondR.ptr(icurr_lapack_rcondR,0),
	     &rcondR,rcondDblWork.ptr(0,0),rcondIntWork.ptr(0,0),
	     &info);
  lapackRcondR(icurr_lapack_rcondR,0)=rcondR;


  int inext_lapack_rcondR=icurr_lapack_rcondR;
  if((rcondR<=min_allowed_rcond)&&
     (inext_lapack_rcondR-iprev_lapack_rcondR==1)) {
    //at this point the previous lapack rcondR==1.0
    rcondR=1.0;
    inext_lapack_rcondR=iprev_lapack_rcondR;
  }

  //do the bisection search if necessary, at most ceil(log2()) more 
  //calls to rcond
  int rcond_iter=0;
  int max_rcond_iter=
    std::ceil(std::log(static_cast<double>
		       (inext_lapack_rcondR-iprev_lapack_rcondR))
	      /std::log(2.0));
  while((lapackRcondR(inext_lapack_rcondR,0)<=min_allowed_rcond)&&
        (inext_lapack_rcondR>iprev_lapack_rcondR)) {
    ++rcond_iter;
    icurr_lapack_rcondR=(iprev_lapack_rcondR+inext_lapack_rcondR)/2;
    numEqnKeep=icurr_lapack_rcondR+1;

    DPOCON_F77(&uplo,&numEqnKeep,RChol.ptr(0,0),&ld_RChol,
	       oneNormPrecondR.ptr(icurr_lapack_rcondR,0),
	       &rcondR,rcondDblWork.ptr(0,0),rcondIntWork.ptr(0,0),
	       &info);
    lapackRcondR(icurr_lapack_rcondR,0)=rcondR;
    //printf("rcond_iter=%d icurr_lapack_rcondR=%d rcondR=%g\n",
    //rcond_iter,icurr_lapack_rcondR,rcondR);

    if(rcondR<min_allowed_rcond)
      inext_lapack_rcondR=icurr_lapack_rcondR;
    else if(min_allowed_rcond<rcondR)
      iprev_lapack_rcondR=icurr_lapack_rcondR;
    else if(min_allowed_rcond==rcondR) {
      //numEqnKeep=icurr_lapack_rcondR+1;
      break;
    }
    if((inext_lapack_rcondR-iprev_lapack_rcondR==1)||
       (max_rcond_iter<rcond_iter)) {
      numEqnKeep=iprev_lapack_rcondR+1;
      rcondR=lapackRcondR(iprev_lapack_rcondR,0);
      break;
    }
  }
#ifdef __TIME_PIVOT_CHOLESKY__
  gettimeofday(&tv, NULL); 
  stop_sec_rcond=tv.tv_sec;
  stop_usec_rcond=tv.tv_usec;
  time_spent_on_rcond_in_pivot_cholesky+=
    (static_cast<double>(stop_sec_rcond -start_sec_rcond )+
     static_cast<double>(stop_usec_rcond-start_usec_rcond)/1000000.0);
#endif

  numRowsR=numEqnKeep;
  polyOrder=polyOrderRequested;
  while((numRowsR<=numTrend(polyOrder,0))&&(polyOrder>0))
    --polyOrder;
  nTrend=numTrend(polyOrder,0);
  Poly.resize(nTrend,numVarsr); //I am relying on the matrix class's actual 
  //size not changing and that it's contents aren't being written over, so 
  //that when I enlarge the matrix up to polyOrderRequested I recover all 
  //the polynomial powers up to polyOrderRequested

  //printf("numRowsR=%d rcondR=%g, time_spent_on_rcond=%g max_rcond_iter=%d rcond_iter=%d\n",
  //numRowsR,rcondR,time_spent_on_rcond_in_pivot_cholesky,max_rcond_iter,rcond_iter);

  RChol.resize(numRowsR,numRowsR); //resize() instead of newSize() 
  //because we want to keep the current contents in the same 2D 
  //order

  Y.newSize(numEqnKeep,1); //newSize() because we don't care about 
  //the current contents of Y
  G.newSize(numEqnKeep,nTrend);
  iEqnKeep.newSize(numEqnKeep,1);
  iptIderKeep.newSize(numEqnKeep,2);
  //ifPointUsed.newSize(numPoints,1);
  //ifPointUsed.zero();
  //int used=1;  

  for(int j=0; j<numEqnKeep; ++j) {
    int jeqn=iEqnKeep(j,0);
    int jpt=jeqn%numPoints;
    //ifPointUsed(jpt,0)=used; 
    iptIderKeep(j,0)=jpt;
    iptIderKeep(j,1)=jeqn/numPoints-1;
    Y(j,0)=Yall(jeqn,0);
  }

  //scale (unprecondition) RChol
  for(int j=0; j<numEqnKeep; ++j) {
    RChol(j,j)*=scaleRChol(iEqnKeep(j,0),1);
    for(int i=j+1; i<numEqnKeep; ++i)
      RChol(j,i)=RChol(i,j)*=scaleRChol(iEqnKeep(i,0),1);
  }

  for(int itrend=0; itrend<nTrend; ++itrend) 
    for(int i=0; i<numEqnKeep; ++i)       
      G(i,itrend)=Gall(iEqnKeep(i,0),itrend);
  

#ifdef __DEBUG_MY_CHOL2__
  {
  FILE* fp=fopen("R_from_eqn_selecting_chol.txt","w");
  /*
  MtxDbl rr;
  correlation_matrix(rr, XR);

  fprintf(fp,"nrows=%d\n\nR=",numEqnKeep);

  for(int i=0; i<numEqnKeep; ++i){
    int iek=iEqnKeep(i,0);
    int jek=iEqnKeep(0,0);
    fprintf(fp,"\n%-22.16g",R(iek,jek)*scaleRChol(iek,1)*scaleRChol(jek,1));
    for(int j=1; j<numEqnKeep; ++j) {
      jek=iEqnKeep(j,0);
      fprintf(fp," %-22.16g",R(iek,jek)*scaleRChol(iek,1)*scaleRChol(jek,1));
    }
  }

  fprintf(fp,"\n\nr=");
  for(int i=0; i<rr.getNRows(); ++i){
    fprintf(fp,"\n%-22.16g",rr(i,0));
    for(int j=1; j<rr.getNCols(); ++j)
      fprintf(fp," %-22.16g",rr(i,j));
  }

  MtxDbl Ichop;
  solve_after_Chol_fact(Ichop,RChol,rr,'T');
  //assert(false);

  fprintf(fp,"\n\nIchop=");  
  for(int i=0; i<Ichop.getNRows(); ++i){
    fprintf(fp,"\n%-22.16g",Ichop(i,0));
    for(int j=1; j<Ichop.getNCols(); ++j)
      fprintf(fp," %-22.16g",Ichop(i,j));
  }

  */

  //RChol.resize(numEqnKeep,numEqnKeep);
  MtxDbl I(numEqnKeep,numEqnKeep);
  MtxDbl RRR(numEqnKeep,numEqnKeep);
  for(int j=0; j<numEqnKeep; ++j) {
    int jek=iEqnKeep(j,0);
    for(int i=0; i<numEqnKeep; ++i) {
      int iek=iEqnKeep(i,0);
      RRR(i,j)=R(iek,jek)*scaleRChol(iek,1)*scaleRChol(jek,1);
    }
  }
  solve_after_Chol_fact(I,RChol,RRR);

  //Rinv.copy(rchol);
  //Rinv.newSize(numEqnKeep,numEqnKeep);
  //inverse_after_Chol_fact(Rinv);

  fprintf(fp,"\n\nRChol=");
  for(int i=0; i<numEqnKeep; ++i) {
    fprintf(fp,"\n");
    for(int j=0; j<=i; ++j) 
      fprintf(fp,"%-12.6g ",RChol(i,j)); 
  }

  fprintf(fp,"\n\nI=");
  for(int i=0; i<numEqnKeep; ++i) {
    fprintf(fp,"\n");
    for(int j=0; j<numEqnKeep; ++j) 
      fprintf(fp,"%-22.16g ",I(i,j));
  }
  /*
  fprintf(fp,"\n\nI=");
  for(int i=0; i<numEqnKeep; ++i) {
    fprintf(fp,"\n");
    for(int j=0; j<numEqnKeep; ++j) {
      int jek=iEqnKeep(j,0);
      tempdouble2=0.0;
      for(int k=0; k<numEqnKeep; ++k) {
	int kek=iEqnKeep(k,0);
	tempdouble2+=Rinv(i,k)*R(kek,jek);
      }
      fprintf(fp,"%-22.16g ",tempdouble2);
    }
  }
  */
  fflush(fp);
  fclose(fp);

  //assert(false);
  }
#endif
    
#ifdef __DEBUG_MY_CHOL__
  {
  FILE* fp=fopen("Check_my_precond_chol2.txt","w");
  fprintf(fp,"nrows=%d\n\nR=",numEqnKeep);

  for(int i=0; i<numEqnKeep; ++i){
    fprintf(fp,"\n%-22.16g",R(i,0));
    for(int j=1; j<numEqnKeep; ++j)
      fprintf(fp," %-22.16g",R(i,j));
  }
  
  double rcondR_ref;
  MtxDbl RCholDebug(R);
  Chol_fact(RCholDebug, info_local, rcondR_ref);

  fprintf(fp,"\n\nrcondR_ref=%-22.16g\nrcondR_my=%-22.16g\n\nL_ref=",
	  rcondR_ref,rcondR);

  for(int i=0; i<numEqnKeep; ++i){
    fprintf(fp,"\n%-22.16g",RCholDebug(i,0));
    for(int j=1; j<=i; ++j)
      fprintf(fp," %-22.16g",RCholDebug(i,j));
    for(int j=i+1; j<numEqnKeep; ++j)
      fprintf(fp," %-22.16g",0.0);
  }  
  
  fprintf(fp,"\n\nL_my=");
  
  for(int i=0; i<numEqnKeep; ++i){
    fprintf(fp,"\n%-22.16g",RChol(i,0));
    for(int j=1; j<=i; ++j)
      fprintf(fp," %-22.16g",RChol(i,j));
    for(int j=i+1; j<numEqnKeep; ++j)
      fprintf(fp," %-22.16g",0.0);
  }  

  fclose(fp);  
  }
  //assert(false);
#endif
#ifdef __TIME_PIVOT_CHOLESKY__
  gettimeofday(&tv, NULL); 
  stop_sec_cholesky=tv.tv_sec;
  stop_usec_cholesky=tv.tv_usec;
  time_spent_on_pivot_cholesky+=
    (static_cast<double>(stop_sec_cholesky -start_sec_cholesky )+
     static_cast<double>(stop_usec_cholesky-start_usec_cholesky)/1000000.0);
  //printf("GEK: numVarsr=%d numPoints=%d numEqnAvail=%d numEqnKeep=%d\ntime_spent_on_pivot_cholesky=%g,rcondR=%g\n",  
  //numVarsr, numPoints, numEqnAvail, numEqnKeep,
  // time_spent_on_pivot_cholesky,rcondR);
  //exit(0);
#endif
  //exit(0);
  return;
}
#elif defined(__PIVOT_CHOL_GEK_R__) 
void GradKrigingModel::equationSelectingPrecondCholR(){ 
#ifdef __TIME_PIVOT_CHOLESKY__
  ++n_pivot_cholesky_calls;
  struct timeval tv;
  gettimeofday(&tv, NULL); 
  long int start_sec_block4, start_usec_block4;
  long int start_sec_cholesky=tv.tv_sec;
  long int start_usec_cholesky=tv.tv_usec;
  long int stop_sec_cholesky, stop_usec_cholesky;
  long int start_sec_rcond, start_usec_rcond;
  long int stop_sec_rcond, stop_usec_rcond;
#endif  
  double min_allowed_rcond=1.0/maxCondNum;
  //double min_allowed_pivot_est_rcond=256.0/maxCondNum;

  scaleRChol.newSize(numEqnAvail,2); 
  for(int i=0; i<numEqnAvail; ++i) {
    scaleRChol(i,1)=std::sqrt(R(i,i));
    scaleRChol(i,0)=1.0/scaleRChol(i,1);
  }
  //precondition R
  for(int j=0; j<numEqnAvail; ++j) {
    for(int i=0; i<numEqnAvail; ++i)
      R(i,j)*=scaleRChol(i,0)*scaleRChol(j,0);
    R(j,j)=1.0; //round off error on the diagonal can screw up pivoting 
    //cholesky... since there is zero correlation between an individual 
    //point's function value and its derivatives and since there are 
    //only ones on the diagonal of the preconditioned R, then pivoting
    //cholesky will not pivot the first point's function value or 
    //derivatives (because there won't be any diagonal entries larger 
    //than one, and the preconditioned correlation between first point's 
    //function value and derivatives is the identity matrix) round off
    //error on the diagonal can screw with that and I rely on those 
    //assumptions to ensure that the anchor point's derivatives won't be
    //"dropped"
  }

  RChol.newSize(numEqnAvail,numEqnAvail);
  int ld_RChol=RChol.getNRowsAct();

  //if the user specifies an anchor point we will do an initial 
  //pivoting cholesky with the anchor point's function value and 
  //derivatives first, followed by the function values only at 
  //the rest of the data points
  //
  //if the user doesn't specify an anchor point we will do an 
  //initial pivoting with only the function values
  //
  //the purpose is to pivot/order the points in a good way, i.e.
  //by pivoting Cholesky, based on a small number of equations, 
  //rather than all equations, because pivoting Cholesky is 
  //significantly more expensive than the highly optimized 
  //non-pivoting LAPACK Cholesky which (if needed) will be 
  //applied to the full set of equations after the points have 
  //been ordered according to the previous pivoting.
  //
  //because the function value and its derivatives at an 
  //individual point are uncorrelated with each other and we
  //have normalized R so it has all ones on it's diagonals,
  //pivoting Cholesky WILL NOT PIVOT the first point's 
  //function values or derivatives (the first equation that 
  //can be pivoting is the second point's function value, read
  //as the firs equation that does not have zero correlation
  //with the first point or any of its derivatives)

  int num_eqn_in_grp=numVarsr+1; 
  int num_eqn_test=0;
  iOrderEqnTest.newSize(numEqnAvail,1);

  if((ifHaveAnchorPoint==true)&&(iAnchorPoint!=0)){
    iOrderEqnTest(0,0)=iAnchorPoint;
    int ieqn=1;
    for(int ivar=0; ivar<numVarsr; ++ivar, ++ieqn) 
      iOrderEqnTest(ieqn,0)=iAnchorPoint+(ivar+1)*numPoints;      
    for(int ipt=0; ipt<iAnchorPoint; ++ipt) {
      iOrderEqnTest(ieqn,0)=ipt;
      ++ieqn;
      for(int ivar=0; ivar<numVarsr; ++ivar, ++ieqn)
	iOrderEqnTest(ieqn,0)=ipt+(ivar+1)*numPoints;
    }
    for(int ipt=iAnchorPoint+1; ipt<numPoints; ++ipt) {
      iOrderEqnTest(ieqn,0)=ipt;
      ++ieqn;
      for(int ivar=0; ivar<numVarsr; ++ivar, ++ieqn)
	iOrderEqnTest(ieqn,0)=ipt+(ivar+1)*numPoints;
    }
  }
  else{
    int ieqn=0;
    for(int ipt=0; ipt<numPoints; ++ipt) {
      iOrderEqnTest(ieqn,0)=ipt;
      ++ieqn;
      for(int ivar=0; ivar<numVarsr; ++ivar, ++ieqn)
	iOrderEqnTest(ieqn,0)=ipt+(ivar+1)*numPoints;
    }
  }

  //copy reordered R to RChol
  for(int j=0; j<numEqnAvail; ++j) {
    int jeqn=iOrderEqnTest(j,0);
    for(int i=0; i<numEqnAvail; ++i) 
      RChol(i,j)=R(iOrderEqnTest(i,0),jeqn);
  }

  int info=0;
  char uplo='B'; //'B' means we have both halves of R in RChol so the 
  //fortran doesn't have to copy one half to the other
  iPivot.newSize(numEqnAvail,1);  
  numEqnKeep=numEqnAvail;
  int blocksize=1+numVarsr;
  MtxDbl dwork(blocksize,3*blocksize);
  BLOCKPIVOTCHOL_F77(&uplo, &numEqnAvail, RChol.ptr(0,0), &ld_RChol,
		     &blocksize, iPivot.ptr(0,0), &numEqnKeep, dwork.ptr(0,0),
		     &min_allowed_rcond, &info); 
  //printf("numEqnAvail=%d numEqnKeep=%d\n",numEqnAvail,numEqnKeep);
  //for(int i=0; i<numEqnKeep; ++i) {
  //printf("i=%d iPivot=%d\n",i,iPivot(i,0));
  //}
  //exit(0);
  
  iEqnKeep.newSize(numEqnAvail,1);
  for(int i=0; i<numEqnKeep; ++i) {
    iEqnKeep(i,0)=iOrderEqnTest(iPivot(i,0)-1,0);
    //printf("i=%d iPivot=%d iEqnKeep=%d\n",i,iPivot(i,0),iEqnKeep(i,0));
  }

  lapackRcondR.newSize(numEqnAvail,1);
  lapackRcondR.zero(); 
  //the normalized correlation matrix for 1 point (function value plus first
  //derivatives) is the identity matrix which has a condition number of 1
  for(int i=0; i<numVarsr+1; ++i)
    lapackRcondR(i,0)=1.0;
  int iprev_lapack_rcondR=numVarsr;

  oneNormPrecondR.newSize(numEqnAvail,1);
  sumAbsColPrecondR.newSize(numEqnAvail,1); //used in computing the one norm 

  //precompute and store the one norm after adding every equation 
  //in the pivoted Cholesky whole-point at a time order (where a 
  //whole-point is the function value and its derivatives at an
  //individual point) this is prepwork for the bisection search for the
  //last equation for which the reordered R is not too ill-conditioned
  for(int i=0; i<numEqnKeep; ++i)
    sumAbsColPrecondR(i,0)=std::fabs(R(iEqnKeep(i,0),0));
  oneNormPrecondR(0,0)=sumAbsColPrecondR(0,0);
  
  double tempdouble;
  for(int j=1; j<numEqnKeep; ++j) {
    int jeqn=iEqnKeep(j,0);
    for(int i=0; i<numEqnKeep; ++i)
      sumAbsColPrecondR(i,0)+=std::fabs((iEqnKeep(i,0),jeqn));
    tempdouble=sumAbsColPrecondR(0,0);
    for(int i=1; i<=j; ++i) 
      if(tempdouble<sumAbsColPrecondR(i,0))
	tempdouble=sumAbsColPrecondR(i,0);
    oneNormPrecondR(j,0)=tempdouble;
  }

  //find the last equation that isn't too poorly conditioned, but 
  //first we need the rcond at the point the Cholesky decomposition 
  //terminated.
  uplo='L';
  rcondDblWork.newSize(3*ld_RChol,1);
  rcondIntWork.newSize(ld_RChol,1);
  int icurr_lapack_rcondR=numEqnKeep-1;
  DPOCON_F77(&uplo,&numEqnKeep,RChol.ptr(0,0),&ld_RChol,
	     oneNormPrecondR.ptr(icurr_lapack_rcondR,0),
	     &rcondR,rcondDblWork.ptr(0,0),rcondIntWork.ptr(0,0),
	     &info);
  lapackRcondR(icurr_lapack_rcondR,0)=rcondR;
  //printf("ilasteqn=%d rcondR=%g\n",icurr_lapack_rcondR,rcondR);

  int inext_lapack_rcondR=icurr_lapack_rcondR;
  if((rcondR<=min_allowed_rcond)&&
     (inext_lapack_rcondR-iprev_lapack_rcondR==1)) {
    //at this point the previous lapack rcondR==1.0
    rcondR=1.0;
    inext_lapack_rcondR=iprev_lapack_rcondR;
  }

  //do the bisection search if necessary, at most ceil(log2()) more 
  //calls to rcond
  int rcond_iter=0;
  int max_rcond_iter=
    std::ceil(std::log(static_cast<double>
		       (inext_lapack_rcondR-iprev_lapack_rcondR))
	      /std::log(2.0));
  while((lapackRcondR(inext_lapack_rcondR,0)<=min_allowed_rcond)&&
        (inext_lapack_rcondR>iprev_lapack_rcondR)) {
    ++rcond_iter;
    icurr_lapack_rcondR=(iprev_lapack_rcondR+inext_lapack_rcondR)/2;
    numEqnKeep=icurr_lapack_rcondR+1;

    DPOCON_F77(&uplo,&numEqnKeep,RChol.ptr(0,0),&ld_RChol,
	       oneNormPrecondR.ptr(icurr_lapack_rcondR,0),
	       &rcondR,rcondDblWork.ptr(0,0),rcondIntWork.ptr(0,0),
	       &info);
    lapackRcondR(icurr_lapack_rcondR,0)=rcondR;
    //printf("rcond_iter=%d icurr_lapack_rcondR=%d rcondR=%g\n",
    //rcond_iter,icurr_lapack_rcondR,rcondR);

    if(rcondR<min_allowed_rcond)
      inext_lapack_rcondR=icurr_lapack_rcondR;
    else if(min_allowed_rcond<rcondR)
      iprev_lapack_rcondR=icurr_lapack_rcondR;
    else if(min_allowed_rcond==rcondR) {
      //numEqnKeep=icurr_lapack_rcondR+1;
      break;
    }
    if((inext_lapack_rcondR-iprev_lapack_rcondR==1)||
       (max_rcond_iter<rcond_iter)) {
      numEqnKeep=iprev_lapack_rcondR+1;
      rcondR=lapackRcondR(iprev_lapack_rcondR,0);
      break;
    }
  }
#ifdef __TIME_PIVOT_CHOLESKY__
  gettimeofday(&tv, NULL); 
  stop_sec_rcond=tv.tv_sec;
  stop_usec_rcond=tv.tv_usec;
  time_spent_on_rcond_in_pivot_cholesky+=
    (static_cast<double>(stop_sec_rcond -start_sec_rcond )+
     static_cast<double>(stop_usec_rcond-start_usec_rcond)/1000000.0);
#endif
  numRowsR=numEqnKeep;
  //printf("numRowsR=%d rcondR=%g, time_spent_on_rcond=%g max_rcond_iter=%d rcond_iter=%d\n",
  //numRowsR,rcondR,time_spent_on_rcond_in_pivot_cholesky,max_rcond_iter,rcond_iter);

  polyOrder=polyOrderRequested;
  while((numRowsR<=numTrend(polyOrder,0))&&(polyOrder>0))
    --polyOrder;
  nTrend=numTrend(polyOrder,0);
  Poly.resize(nTrend,numVarsr); //I am relying on the matrix class's actual 
  //size not changing and that it's contents aren't being written over, so 
  //that when I enlarge the matrix up to polyOrderRequested I recover all 
  //the polynomial powers up to polyOrderRequested

  RChol.resize(numRowsR,numRowsR); //resize() instead of newSize() 
  //because we want to keep the current contents in the same 2D 
  //order

  Y.newSize(numEqnKeep,1); //newSize() because we don't care about 
  //the current contents of Y
  //int ntrend=getNTrend();
  G.newSize(numEqnKeep,nTrend);
  iEqnKeep.newSize(numEqnKeep,1);
  iptIderKeep.newSize(numEqnKeep,2);
  //ifPointUsed.newSize(numPoints,1);
  //ifPointUsed.zero();
  //int used=1;  

  for(int j=0; j<numEqnKeep; ++j) {
    int jeqn=iEqnKeep(j,0);
    int jpt=jeqn%numPoints;
    //ifPointUsed(jpt,0)=used; 
    iptIderKeep(j,0)=jpt;
    iptIderKeep(j,1)=jeqn/numPoints-1;
    Y(j,0)=Yall(jeqn,0);
  }

  //scale (unprecondition) RChol
  for(int j=0; j<numEqnKeep; ++j) 
    for(int i=j; i<numEqnKeep; ++i)
      RChol(i,j)*=scaleRChol(iEqnKeep(i,0),1);
  

  for(int itrend=0; itrend<nTrend; ++itrend) 
    for(int i=0; i<numEqnKeep; ++i) {      
      int ieqn=iEqnKeep(i,0);
      G(i,itrend)=Gall(ieqn,itrend);
    }
  

#ifdef __DEBUG_MY_CHOL2__
  {
  FILE* fp=fopen("R_from_eqn_selecting_chol.txt","w");
  fprintf(fp,"\n\nRChol=");
  for(int i=0; i<numEqnKeep; ++i) {
    fprintf(fp,"\n%-2d ",i);
    for(int j=0; j<=i; ++j) 
      fprintf(fp,"%-12.6g ",RChol(i,j)); 
  }

  /*
  MtxDbl rr;
  correlation_matrix(rr, XR);

  fprintf(fp,"nrows=%d\n\nR=",numEqnKeep);

  for(int i=0; i<numEqnKeep; ++i){
    int iek=iEqnKeep(i,0);
    int jek=iEqnKeep(0,0);
    fprintf(fp,"\n%-22.16g",R(iek,jek)*scaleRChol(iek,1)*scaleRChol(jek,1));
    for(int j=1; j<numEqnKeep; ++j) {
      jek=iEqnKeep(j,0);
      fprintf(fp," %-22.16g",R(iek,jek)*scaleRChol(iek,1)*scaleRChol(jek,1));
    }
  }

  fprintf(fp,"\n\nr=");
  for(int i=0; i<rr.getNRows(); ++i){
    fprintf(fp,"\n%-22.16g",rr(i,0));
    for(int j=1; j<rr.getNCols(); ++j)
      fprintf(fp," %-22.16g",rr(i,j));
  }

  MtxDbl Ichop;
  solve_after_Chol_fact(Ichop,RChol,rr,'T');
  //assert(false);

  fprintf(fp,"\n\nIchop=");  
  for(int i=0; i<Ichop.getNRows(); ++i){
    fprintf(fp,"\n%-22.16g",Ichop(i,0));
    for(int j=1; j<Ichop.getNCols(); ++j)
      fprintf(fp," %-22.16g",Ichop(i,j));
  }

  */

  //RChol.resize(numEqnKeep,numEqnKeep);
  MtxDbl I(numEqnKeep,numEqnKeep);
  MtxDbl RRR(numEqnKeep,numEqnKeep);
  for(int j=0; j<numEqnKeep; ++j) {
    int jek=iEqnKeep(j,0);
    for(int i=0; i<numEqnKeep; ++i) {
      int iek=iEqnKeep(i,0);
      RRR(i,j)=R(iek,jek)*scaleRChol(iek,1)*scaleRChol(jek,1);
    }
  }
  int inforrr;
  double rcondRRR;
  MtxDbl rrrchol(RRR);
  Chol_fact(rrrchol,inforrr,rcondRRR);
  solve_after_Chol_fact(I,RChol,RRR);

  //Rinv.copy(rchol);
  //Rinv.newSize(numEqnKeep,numEqnKeep);
  //inverse_after_Chol_fact(Rinv);

  fprintf(fp,"\n\nrrrChol=");
  for(int i=0; i<numEqnKeep; ++i) {
    fprintf(fp,"\n%-2d ",i);
    for(int j=0; j<=i; ++j) 
      fprintf(fp,"%-12.6g ",rrrchol(i,j)); 
  }

  for(int j=0; j<numEqnKeep; ++j)
    for(int i=0; i<numEqnKeep; ++i)
      R(i,j)*=scaleRChol(i,1)*scaleRChol(j,1);
  Chol_fact(R,inforrr,rcondRRR);
  
  fprintf(fp,"\n\nChol(R)=");
  for(int i=0; i<numEqnKeep; ++i) {
    fprintf(fp,"\n%-2d ",i);
    for(int j=0; j<=i; ++j) 
      fprintf(fp,"%-12.6g ",R(i,j)); 
  }

  fprintf(fp,"\n\nI=");
  for(int i=0; i<numEqnKeep; ++i) {
    fprintf(fp,"\n%-2d ",i);
    for(int j=0; j<numEqnKeep; ++j) 
      fprintf(fp,"%-22.16g ",I(i,j));
  }
  /*
  fprintf(fp,"\n\nI=");
  for(int i=0; i<numEqnKeep; ++i) {
    fprintf(fp,"\n");
    for(int j=0; j<numEqnKeep; ++j) {
      int jek=iEqnKeep(j,0);
      tempdouble2=0.0;
      for(int k=0; k<numEqnKeep; ++k) {
	int kek=iEqnKeep(k,0);
	tempdouble2+=Rinv(i,k)*R(kek,jek);
      }
      fprintf(fp,"%-22.16g ",tempdouble2);
    }
  }
  */
  fflush(fp);
  fclose(fp);

  assert(false);
  }
#endif
    
#ifdef __DEBUG_MY_CHOL__
  {
  FILE* fp=fopen("Check_my_precond_chol2.txt","w");
  fprintf(fp,"nrows=%d\n\nR=",numEqnKeep);

  for(int i=0; i<numEqnKeep; ++i){
    fprintf(fp,"\n%-22.16g",R(i,0));
    for(int j=1; j<numEqnKeep; ++j)
      fprintf(fp," %-22.16g",R(i,j));
  }
  
  double rcondR_ref;
  MtxDbl RCholDebug(R);
  Chol_fact(RCholDebug, info_local, rcondR_ref);

  fprintf(fp,"\n\nrcondR_ref=%-22.16g\nrcondR_my=%-22.16g\n\nL_ref=",
	  rcondR_ref,rcondR);

  for(int i=0; i<numEqnKeep; ++i){
    fprintf(fp,"\n%-22.16g",RCholDebug(i,0));
    for(int j=1; j<=i; ++j)
      fprintf(fp," %-22.16g",RCholDebug(i,j));
    for(int j=i+1; j<numEqnKeep; ++j)
      fprintf(fp," %-22.16g",0.0);
  }  
  
  fprintf(fp,"\n\nL_my=");
  
  for(int i=0; i<numEqnKeep; ++i){
    fprintf(fp,"\n%-22.16g",RChol(i,0));
    for(int j=1; j<=i; ++j)
      fprintf(fp," %-22.16g",RChol(i,j));
    for(int j=i+1; j<numEqnKeep; ++j)
      fprintf(fp," %-22.16g",0.0);
  }  

  fclose(fp);  
  }
  //assert(false);
#endif
#ifdef __TIME_PIVOT_CHOLESKY__
  gettimeofday(&tv, NULL); 
  stop_sec_cholesky=tv.tv_sec;
  stop_usec_cholesky=tv.tv_usec;
  time_spent_on_pivot_cholesky+=
    (static_cast<double>(stop_sec_cholesky -start_sec_cholesky )+
     static_cast<double>(stop_usec_cholesky-start_usec_cholesky)/1000000.0);
  //printf("GEK: numVarsr=%d numPoints=%d numEqnAvail=%d numEqnKeep=%d\ntime_spent_on_pivot_cholesky=%g,rcondR=%g\n",  
  //numVarsr, numPoints, numEqnAvail, numEqnKeep,
  // time_spent_on_pivot_cholesky,rcondR);
  //exit(0);
#endif
  //exit(0);
  return;
}
#endif




/***********************************************************************/
/***********************************************************************/
/**** Unit Test functions for Kriging start here                    ****/
/***********************************************************************/
/***********************************************************************/




/***********************************************************************/
/***********************************************************************/
/**** Unit Test functions for Kriging end here                      ****/
/***********************************************************************/
/***********************************************************************/


/***********************************************************************/
/***********************************************************************/
/**** KrigingModelFactory member functions start here               ****/
/***********************************************************************/
/***********************************************************************/

/***********************************************************************/
/***********************************************************************/
/**** KrigingModelFactory member functions end here                 ****/
/***********************************************************************/
/***********************************************************************/

/***********************************************************************/
/***********************************************************************/
/**** KrigingModel member functions start here                      ****/
/***********************************************************************/
/***********************************************************************/

// typical constructor
GradKrigingModel::GradKrigingModel(const SurfData& sd, const ParamMap& params)
  : SurfPackModel(sd,sd.getJOut()), numVarsr(sd.getNVarsr()), 
    numTheta(numVarsr), numPoints(sd.getNPts()), XR(sdBuild.xr)
{
#ifdef __TIME_PIVOT_CHOLESKY__
  n_pivot_cholesky_calls=0;
  n_rcond_calls_in_pivot_cholesky=0;
  time_spent_on_rcond_in_pivot_cholesky=0.0;
  time_spent_on_pivot_cholesky=0.0;
  time_spent_on_pivot_cholesky_block1=0.0;
  time_spent_on_pivot_cholesky_blocks1_2=0.0;
  time_spent_on_pivot_cholesky_blocks1_2_3=0.0;
  time_spent_on_pivot_cholesky_block4=0.0;
#endif

  multi_dim_poly_power(Der, numVarsr, 1);  //use all mixed partial 
  //derivatives, up to first order, of the basis functions
  nDer=Der.getNRows(); //for GradKrigingModel nDer=(1+numVarsr);

  //numRowsR=numPoints*nDer;
  //printf("calling the right GradKrigingModel constructor\n"); fflush(stdout);
  
  numEqnAvail=numPoints*nDer;

  //if the SurfDataScaler class does what it's supposed to (the only private content in sdBuild that a Model can access are the scaling bits, and only the model can see inside the scaler) the next line will cause an error when you try to compile with it uncommented
  //printf("sd.scaler->mySd.jout=%d\n",scaler.mySd.jout);
  
  // OPTIONS PARSING
  // BMA TODO: max_trials, lower/upper bounds (use toVec helper),
  //           correlation_lengths, optimizationMethod
  ParamMap::const_iterator param_it;

  //printf("GradKrigingModel() says hello world\n");


  // *************************************************************  
  // control verbosity outputLevel
  // *************************************************************  
  param_it = params.find("verbosity");
  if (param_it != params.end() && param_it->second.size() > 0)
    outputLevel=static_cast<short>(std::atoi(param_it->second.c_str()));
  
  // *************************************************************  
  // detect an anchor point if present this is the one point that 
  // we make sure that the equationSelectingCholR does not discard
  // *************************************************************
  iAnchorPoint=-1;
  ifHaveAnchorPoint=false;
  //std::cerr << "A: ifHaveAnchorPoint=" << ifHaveAnchorPoint <<
  //" iAnchorPoint=" << iAnchorPoint << std::endl;
  param_it = params.find("anchor_index");
  if (param_it != params.end() && param_it->second.size() > 0) {
    ifHaveAnchorPoint=true;
    //printf("GradKrigingModel() sees an anchor point\n");
    //fflush(stdout);
    iAnchorPoint=std::atoi(param_it->second.c_str());
    //printf("iAnchorPoint=%d\n",iAnchorPoint);
    //fflush(stdout);
    assert((0<=iAnchorPoint)&&(iAnchorPoint<numPoints));
  }
  //std::cerr << "B: ifHaveAnchorPoint=" << ifHaveAnchorPoint <<
  //" iAnchorPoint=" << iAnchorPoint << std::endl;
  //fflush(stderr);

  //may rename or remove this, modify for a single anchor point
  getBaseIEqnKeep();



  //assert(false)
  // *************************************************************
  // this starts the input section about scaling the data
  // *************************************************************

  MtxDbl min_max_xr(2,numVarsr);
  bool if_user_specified_lower_bounds=false;
  param_it = params.find("lower_bounds");
  if (param_it != params.end() && param_it->second.size() > 0) {
    if_user_specified_lower_bounds=true;
    assert(!(min_max_xr.putRows(param_it->second,0)));
  }
  
  bool if_user_specified_upper_bounds=false;
  param_it = params.find("upper_bounds");
  if (param_it != params.end() && param_it->second.size() > 0) {
    if_user_specified_upper_bounds=true;
    assert(!(min_max_xr.putRows(param_it->second,1)));
  }
  
  if(!(if_user_specified_lower_bounds==if_user_specified_upper_bounds)) {
    std::cerr << "Your options are to\n(A) specify both the upper and lower, or\n(B) specify neither the upper nor lower,\nbounds of the domain of the Kriging Model\n";
    assert(if_user_specified_lower_bounds==if_user_specified_upper_bounds);
  }
  
  if(if_user_specified_lower_bounds==true) {
    for(int ivarsr=0; ivarsr<numVarsr; ++ivarsr) 
      if(!(min_max_xr(0,ivarsr)<=min_max_xr(1,ivarsr))) {
	std::cerr << "The lower bound of the domain of the Kriging Model must be less than or equal to the upper bound of the domain of the Kriging Model\n";
	assert(min_max_xr(0,ivarsr)<=min_max_xr(1,ivarsr));
      }
    //printf("lower_bounds = (%g",min_max_xr(0,0));
    //for(int ivarsr=1; ivarsr<numVarsr; ++ivarsr)
    //printf(", %g",min_max_xr(0,ivarsr));
    //printf("), upper_bounds = (%g",min_max_xr(1,0));
    //for(int ivarsr=1; ivarsr<numVarsr; ++ivarsr)
    //printf(", %g",min_max_xr(1,ivarsr));    
    //printf(")\n");
    sdBuild.setUnscaledDomainSize(min_max_xr);
  }
  
  //printf("GradKrigingModel constructor should have just written out domain bounds\n");
  
  param_it = params.find("dimension_groups");
  if (param_it != params.end() && param_it->second.size() > 0) {
    MtxInt dim_groups(1,numVarsr);
    assert(!(dim_groups.putRows(param_it->second,0)));
    sdBuild.setDimGroups(dim_groups);
  }
  
  scaler.scaleToDefault(); //scale outputs to -0.5<=Y<=0.5 and scale real inputs to volume 1 hyper-rectangle centered at 0 if real iputs dimensions are locked or the unit hypercube centered at 0 if no dimensions are locked.  The scaling is done to let us define the feasible region simply (done in create);
  
  //std::string debug_filename = "GradRos10sdBuildScaled.spd";
  //sdBuild.write(debug_filename);


  sdBuild.getUpToDerY(Yall,1); //Y contains [ sdBuild.y(:,jout) sdBuild.derY[jout][1] ] (the notation of this comment mixes MATLAB commands with SurfData variables)
  Yall.reshape(numEqnAvail,1);


  //printf("size(Y)=[%d %d]\n",Y.getNRows(),Y.getNCols());


  // *************************************************************
  // this starts the input section about optimizing or directly
  // scepcifying correlation lengths, it must come after the 
  // scaling section
  // *************************************************************
  
  // current options are none (fixed correl) | sampling (guess) | local | global | global_local
  optimizationMethod = "global_local";
  param_it = params.find("optimization_method");
  if (param_it != params.end() && param_it->second.size() > 0)
    optimizationMethod = param_it->second; 
  
  if(optimizationMethod.compare("none")==0)
    maxTrials=1;
  else if(optimizationMethod.compare("local")==0)
    maxTrials=20;
  else if(optimizationMethod.compare("sampling")==0)
    maxTrials=2*numVarsr+1;
  else if(optimizationMethod.compare("global")==0) {
    //maxTrials=(2*numVarsr+1)*(numVarsr+1)*10;
    //if(maxTrials>1500) 
    //maxTrials=1500;
    maxTrials = 10000;
  }
  else if(optimizationMethod.compare("global_local")==0) {
    maxTrials = 10000; //ensure it has non-zero as a fail safe but this shouldn't be used
    maxTrialsGlobal = 500;
    maxTrialsLocal = 20;
  }
  else{ //error checking the input
    std::cerr << "GradKrigingModel() unknown optimization_method [" << optimizationMethod << "]  aborting\n";
    assert(false);
  }

  //cout << "optimization_method=\"" << optimizationMethod << "\"\n";

  //numStarts is the number of starting locations in a multi-start local search
  numStarts=1;
  param_it = params.find("num_starts");
  if (param_it != params.end() && param_it->second.size() > 0) {
    numStarts = std::atoi(param_it->second.c_str()); 
    assert(numStarts>=1);
  }
  
  if(!((numStarts==1)||(optimizationMethod.compare("local")==0))) {
    std::cerr << "Local optimization is the only optimization method for Kriging that uses the \"num_starts\" key word. Check your input file for errors.\n";
    assert((numStarts==1)||(optimizationMethod.compare("local")==0));
  }
  
  //cout << "num_starts=" << numStarts << "\n";


  // does the user want to specify correlation lengths directly?
  ifUserSpecifiedCorrLengths=false; //the default is no
  param_it = params.find("correlation_lengths");
  if (param_it != params.end() && param_it->second.size() > 0) {
    ifUserSpecifiedCorrLengths=true;
    //printf("User specifying correlation lengths\n"); fflush(stdout);
    
    // make sure that the user didn't 
    // * say they want to global optimize __AND__
    // * specify correlation lengths  
    if(optimizationMethod.compare("global")==0) {
      std::cerr << "You can't both \n (A) use the global optimization method to choose, and \n (B) directly specify \n correlation lengths for the Kriging model.\n";
      assert(optimizationMethod.compare("global")!=0);
    }
    else if(optimizationMethod.compare("sampling")==0) {
      // this is only the default number of samples/maxTrials, the user can 
      // still overide this below
      maxTrials+=1;
    }
    
    natLogCorrLen.newSize(1,numVarsr); //allocate space 

    //printf("numVarsr=%d\n",numVarsr); fflush(stdout);
    
    //read the correlation lengths in from the string
    assert(!(natLogCorrLen.putRows(param_it->second,0)));
    // "natLogCorrLen" currently holds the unscaled correlation LENGTHS, not 
    // the natural log of the scaled correlation length, we need to fix that
    // but first we need to check the input for errors
    for(int ivarsr=0; ivarsr<numVarsr; ++ivarsr) 
      if(!(natLogCorrLen(0,ivarsr)>0.0)) {
	std::cerr << "For the Kriging Model, correlation lengths must be strictly positive\n.";
	assert(false);
      }

    //printf("unscaled corr lens = [%12.6g",natLogCorrLen(0,0)); 
    //for(int ivarsr=1; ivarsr<numVarsr; ++ivarsr)
    //printf(", %12.6g",natLogCorrLen(0,ivarsr));
    //printf("]\n");    

    scaler.scaleXrDist(natLogCorrLen); //scale the lengths
    //scaler.scaleXrOther(natLogCorrLen); //error
    //printf("scaled corr lens = [%12.6g",natLogCorrLen(0,0)); 
    //for(int ivarsr=1; ivarsr<numVarsr; ++ivarsr)
    // printf(", %12.6g",natLogCorrLen(0,ivarsr));
    //printf("]\n");    
    //fflush(stdout);
    
    //compute the natural log of the correlation lengths
    for(int ivarsr=0; ivarsr<numVarsr; ++ivarsr) 
      natLogCorrLen(0,ivarsr)=std::log(natLogCorrLen(0,ivarsr)); 
    
    //natLogCorrLen will be the first of the initial iterates (guesses), this happens in the create() function below
  }
  //printf("If user specified correlationlengths we should have just printed them\n");

  // maximum objective evals for optimization or guess
  param_it = params.find("max_trials");
  if (param_it != params.end() && param_it->second.size() > 0) {
    maxTrials = std::atoi(param_it->second.c_str()); 
  }

  assert (maxTrials > 0);

  //printf("maxTrials=%d\n",maxTrials);

  
  // *************************************************************
  // this starts the input section about the trend function 
  // *************************************************************
  polyOrderRequested = 2;
  ifReducedPoly=false;
  param_it = params.find("order");
  if (param_it != params.end() && param_it->second.size() > 0) {
    polyOrderRequested = std::atoi(param_it->second.c_str()); 
    assert (polyOrderRequested >= 0);
  }
  else{
    ifReducedPoly=true;
  }
  numTrend.newSize(polyOrderRequested+1,1);

  //polyOrder = 2; //for debug
  //main_effects_poly_power(Poly, numVarsr, polyOrder); //for debug
  //commented out for debug
  param_it = params.find("reduced_polynomial");
  if (param_it != params.end() && param_it->second.size() > 0) 
    if((std::atoi(param_it->second.c_str()))!=0) 
      ifReducedPoly=true;
      
  //cout << "ifReducedPoly=" << ifReducedPoly << "\n";

  if(ifReducedPoly) {
    main_effects_poly_power(Poly, numVarsr, polyOrderRequested);
    for(polyOrder=0; polyOrder<=polyOrderRequested; ++polyOrder) 
      numTrend(polyOrder,0)=polyOrder*numVarsr+1;
  }
  else{
    multi_dim_poly_power(Poly, numVarsr, polyOrderRequested);  
    for(polyOrder=0; polyOrder<=polyOrderRequested; ++polyOrder) 
      numTrend(polyOrder,0)=num_multi_dim_poly_coef(numVarsr, polyOrder);
  }

  preAllocateMaxMemory(); 
  // *************************************************************
  // this starts the input section HOW to bound the condition 
  // number, this determines which derivatives of the constraint
  // function can be computed analytically so handle that here too
  // *************************************************************
  constraintType="rcond";


  //default orders of analytical derivatives of objective and constraint 
  //functions, eventually (Surfpack release) we will expose this to the
  //user directly
  int num_analytic_obj_ders_in = 0, num_analytic_con_ders_in = 0; 
  //0 = analytic function values
  //1 = analytic function values + first derivatives
  //2 = values + 1st ders + 2nd ders ...
  //currently masterObjectiveAndConstraints only does analytic derivatives
  //up to order 2 (values, and first and second derivatives) of the objective 
  //function. And the method of bounding the condition number determines 
  //the order of analytical derivatives available for the constraints: 
  //eigen values = 1st order analytic derivatives
  //rcond=0th order (values only) analytical derivatives

  //can't do any analytical derivatives (only function values) of rcond
  assert(num_analytic_con_ders_in==0); 
  numConFunc=1;
  
  //convert to the Dakota bitflag convention for derivative orders
  maxObjDerMode=(static_cast<int>(std::pow(2.0,num_analytic_obj_ders_in+1)))-1; //analytical gradients of objective function
  maxConDerMode=(static_cast<int> (std::pow(2.0,num_analytic_con_ders_in+1)))-1; //analytical gradients of constraint function(s)

  //maxCondNum=std::pow(1024.0,5); 
  maxCondNum=std::pow(1024.0,4); 
  //maxCondNum=std::pow(1024.0,5)/32.0;

  // *************************************************************
  // this starts the input section about the nugget which can be
  // used to smooth the data and also decrease the condition 
  // number
  // *************************************************************

  nug=(2*nTrend+1.0)/maxCondNum;
  //nug=2*numPoints/maxCondNum;


  ifChooseNug = false;
  param_it = params.find("find_nugget");
  if (param_it != params.end() && param_it->second.size() > 0)
    ifChooseNug = true; 

  //cout << "ifChooseNug=" << ifChooseNug << "\n";

  // fixed value for now
  maxChooseNug = 0.1;
  nug = 0.0; //default

  nuggetFormula=0;
  param_it = params.find("nugget_formula");
  if (param_it != params.end() && param_it->second.size() > 0) {
    if(ifChooseNug==true) {
      std::cerr << "You can't both auto-select a nugget and use a preset formula" << std::endl;
      assert(ifChooseNug==false);
    }
    nuggetFormula=std::atoi(param_it->second.c_str()); 
    if(nuggetFormula!=0) {
      switch(nuggetFormula) {
      case 1:
	nug=(2*nTrend+1.0)/maxCondNum;
	break;
      case 2:
	nug=2*numPoints/maxCondNum; //bob; may want to change this numPoints to numRowsr
	break;
      default:
	std::cerr << "nugget_formula =" << nuggetFormula << " is not one of the available preset nugget formulas." << std::endl;
	assert(false);
      }
    }
  }

  param_it = params.find("nugget");
  if (param_it != params.end() && param_it->second.size() > 0) {
    if(!((nuggetFormula==0)&&(ifChooseNug==false))) {
      std::cerr << "You can do at most 1 of the following (A) auto-select the nugget (minimum needed to satisfy the condition number bound) (B) use one of the preset nugget formulas (C) directly specify a nugget.  The default is not to use a nugget at all (i.e. use a nugget of zero)." << std::endl;
      assert((nuggetFormula==0)&&(ifChooseNug==false));
    }
    nug = std::atof(param_it->second.c_str()); 
    if(!(nug >= 0.0)) {
      std::cerr << "The nugget must be greater than or equal to zero." << std::endl;
      assert (nug >= 0.0);
    }
  }

  // *************************************************************
  // this ends the input parsing now finish up the prep work
  // *************************************************************

  // BMA TODO: allow user (or at least developer) to specify whether
  // to use rotations or just the identity

  // initialize trend function (and rotations)
  // this initializes EulAng, Rot, Poly, G (trend), and Z (diff on XR)

  //LinearRegressionProblem lrp(sd,polyOrder);
  //poly=lrp.poly;
  //lrp.getInitGuess(EulAng);
  //lrp.optimize_by_picking_best_guess(EulAng); EulAng.reshape(nchoosek(numVarsr, 2),1);
  EulAng.newSize(nchoosek(numVarsr, 2),1); 
  EulAng.zero();
  gen_rot_mat(Rot, EulAng, numVarsr);
  eval_trend_fn(Gall, Poly, Der, Rot, XR);


  /*
  FILE* fp=fopen("GradRos10GMat.txt","w");
  int ntrend=getNTrend();
  for(int i=0; i<numRowsR; ++i) {
    fprintf(fp,"%-12.6g",G(i,0));
    for(int j=1; j<ntrend; ++j)
      fprintf(fp," %-12.6g",G(i,j));
    fprintf(fp,"\n");
  }
  fclose(fp);
  */
  //  assert(false);
    

  //LinearRegressionModel::evalBasis(G,poly,Rot,XR);

  gen_Z_matrix();  

  //printf("completed the right GradKrigingModel constructor\n", stdout); fflush(stdout);
}

void GradKrigingModel::create()
{
  //printf("entered create()\n"); fflush(stdout);

  prevObjDerMode=prevConDerMode=0; //tells us not to reuse previous work used
  //to calculate the objective, constraints and their derivatives the first 
  //time they are asked for
  prevTheta.newSize(1,numTheta); 
  prevTheta.zero(); //not necessary just useful to debug

  //printf("KM.create():1: nug=%g\n",nug);

  // -
  // solve the optimization problem for the correlations
  // -

  //printf("numVarsr=%d\n",numVarsr); fflush(stdout);
  OptimizationProblem opt(*this, numVarsr, numConFunc);

  // set the bounds for the plausible region for correlation lengths
  // (assumes input space has a volume of 1, and data points are
  // uniformly distributed)

  aveDistBetweenPts=std::pow(numPoints,-1.0/numVarsr);
  //aveDistBetweenPts=std::pow(numRowsR,-1.0/numVarsr); //count each derivative as 
  //another point when determining the "average distance between points"
  //this is used to determine the "reach" of each data point...
  //on second thought should I shorten the correlation lengths because I 
  //added derivative information, it might help with conditioning but 
  //should the data points have any less reach because of it?  need to 
  //think about this
  //bob  "bob" is how I flag a section of code that I want to find easily

  /* Gaussian Process error model has about ~5% confidence (2 std devs away) 
     in what points 16 neighbors away have to say. If points are correlated 
     well at even greater distances then either 
     * that same information will be contained in nearby points OR
     * you shouldn't be using a Gaussian correlation function
     KRD */
  double max_corr_length = aveDistBetweenPts*8.0; 

  maxNatLogCorrLen=std::log(max_corr_length);

  /* Gaussian Process error model has about ~5% confidence (2 std devs) midway
     between neighboring points... i.e. you're 4 std devs away from your 
     nearest neighbor so all sample points are treated as being essentially 
     uncorrelated 
     KRD */
  double min_corr_length = aveDistBetweenPts/4.0; 
  minNatLogCorrLen=std::log(min_corr_length);
  //double max_correlation = 1.0/(2.0*min_corr_length*min_corr_length);

  //Choose dead center (in log(correlation length)) of feasible region as the 
  //default initial guess for the Gaussian Process error model, KRD  
  double init_guess=0.5*(maxNatLogCorrLen+minNatLogCorrLen);

  //order the derivative equations for later one pass acceptance/rejectance
  //according to how much new information they provide when correlation lengths
  //at the most favorable (smallest corner of the small feasible region) and 
  //screen out any equations that make it singular to save cholesky decomposition 
  //cost later
  
  MtxDbl theta(1,numTheta);
  theta(0,0)=0.5/(min_corr_length*min_corr_length);
  for(int i=1; i<numTheta; ++i)
    theta(0,i)=theta(0,0);

  //init_guess=minNatLogCorrLen;

  //printf("got to yada yada\n"); fflush(stdout);
  ///set the bounds and the initial iterates
  if(ifUserSpecifiedCorrLengths==true) {
    //printf("says that the user specified correlation lengths\n"); fflush(stdout);
    // the first guess is what the user told us he/she wanted to use
    for(int jvar=0; jvar<numVarsr; ++jvar) {
      opt.lower_bound(jvar, minNatLogCorrLen);
      opt.upper_bound(jvar, maxNatLogCorrLen);
      //double temp_double=natLogCorrLen(0,jvar);
      //printf("GradKrigingModel::create() jvar=%d temp_double=%g\n",jvar,temp_double);
      //fflush(stdout);
      opt.initial_iterate(jvar, natLogCorrLen(0,jvar));
    }
    // the second guess is the center of the small feasible region
    MtxDbl the_second_guess(1,numVarsr);
    for(int jvar=0; jvar<numVarsr; ++jvar) 
      the_second_guess(0,jvar)=init_guess;
    opt.add_initial_iterates(the_second_guess);
  } else {
    
    //printf("about to set bounds and initial iterate\n"); fflush(stdout);
    // since the user didn't specify an initial guess we will use the center
    // of the small feasible region as our first initial guess
    for(int jvar=0; jvar< numVarsr; ++jvar) {
      opt.lower_bound(jvar, minNatLogCorrLen);
      opt.upper_bound(jvar, maxNatLogCorrLen);
      opt.initial_iterate(jvar, init_guess);
    }
  }

  //printf("just set bounds and initial iterate\n"); fflush(stdout);


  //add a bin opt random design with 2*numVars more guesses, 
  //bins are the endpoints of a randomly rotated orthogonal axes
  MtxDbl axes_of_guesses(2*numVarsr,numVarsr);
  gen_rand_axis_bin_opt_samples_0to1(axes_of_guesses,numVarsr);
  for(int i=0; i<2*numVarsr; ++i) {
    //printf("i=%d {",i);
    for(int j=0; j<numVarsr; ++j) {
      //printf(" %g,",axes_of_guesses(i,j));
      axes_of_guesses(i,j)=(maxNatLogCorrLen-minNatLogCorrLen)*axes_of_guesses(i,j)+minNatLogCorrLen;
    }
    //printf("}\n");
  }
  opt.add_initial_iterates(axes_of_guesses);

  //choose the optimizer you want to use
  if(optimizationMethod.compare("none")==0) {
    natLogCorrLen.resize(1,numVarsr);
    opt.retrieve_initial_iterate(0,natLogCorrLen);
  } 
  else{
    if(optimizationMethod.compare("local")==0) {
      if(numStarts==1)
	opt.conmin_optimize();
      else{
	//printf("doing multi-start local optimization\n");
	opt.multistart_conmin_optimize(numStarts);
      }
    }
    else if(optimizationMethod.compare("global")==0)
      opt.direct_optimize();
    else if(optimizationMethod.compare("sampling")==0)
      opt.best_guess_optimize(maxTrials);
    else if(optimizationMethod.compare("global_local")==0){
      maxTrials=maxTrialsGlobal;
      opt.direct_optimize();
      natLogCorrLen = opt.best_point();
      maxTrials=maxTrialsLocal;
      opt.conmin_optimize();
    }
    else{
      std::cerr << "GradKrigingModel:create() unknown optimization_method [" << optimizationMethod << "]  aborting\n";
      assert(false);
    }
    natLogCorrLen = opt.best_point();
  }

  correlations.newSize(1,numVarsr);
  correlations(0,0)=0.5*std::exp(-2.0*natLogCorrLen(0,0));
  //printf("theta={%g",correlations(0,0));
  for(int k=1; k<numVarsr; ++k) {
    correlations(0,k)=0.5*std::exp(-2.0*natLogCorrLen(0,k));
    //printf(", %g",correlations(0,k));
  }

  //printf("}\n");

  //printf("scaled correlations=[%12.6g",correlations(0,0));
  //for(int ivarsr=1; ivarsr<numVarsr; ++ivarsr)
  //printf(", %12.6g",correlations(0,ivarsr));
  //printf("]\n");

  masterObjectiveAndConstraints(correlations, 1, 0);
  if(outputLevel >= NORMAL_OUTPUT)
    std::cout << model_summary_string();
  //deallocate matrices we no longer need after emulator has been created

  //temporary variables used by masterObjectiveAndConstraints
  temp.clear(); //vector
  temp2.clear(); //vector


  //Gtran_Rinv_G_inv.clear(); //need this for derivatives of log(det(Gtran_Rinv_G)) but could use it to replace the permanent copy of Gtran_Rinv_G_Chol

  //variables whose values needed to be retained between sequential call to masterObjectiveAndConstraints for precompute and store strategy to work
  prevObjDerMode=prevConDerMode=0;
  prevTheta.clear(); //row vector 
  Z.clear(); //matrix
  R.clear(); //matrix
  G.clear(); //matrix
  Yall.clear();
  Gall.clear();
  rcondDblWork.clear();
  rcondIntWork.clear();
  sumAbsColPrecondR.clear();
  lapackRcondR.clear();
  oneNormPrecondR.clear();

  con.clear(); //vector
  gradObj.clear(); //vector
  gradCon.clear(); //matrix
  hessObj.clear(); //matrix used to compute hessObj
}

void GradKrigingModel::preAllocateMaxMemory() {
  //this preallocates the maximum size of arrays whose size depends on how many equations were discarded by pivoted Cholesky and they could possibly be allocated to a different size than their maximum the first time they are allocated.
  
  nTrend=numTrend(polyOrderRequested,0);
  Y.newSize(numEqnAvail,1);
  G.newSize(numEqnAvail,nTrend);
  Rinv_G.newSize(numEqnAvail,nTrend);
  Gtran_Rinv_G_Chol.newSize(nTrend,nTrend);
  rhs.newSize(numEqnAvail,1);
  betaHat.newSize(nTrend,1);
  temp.newSize(nTrend,1);
  temp2.newSize(numEqnAvail,1);

  return;
}

std::string GradKrigingModel::model_summary_string() const {
  MtxDbl temp_out_corr_lengths(1,numVarsr);
  for(int i=0; i<numVarsr; ++i) 
    temp_out_corr_lengths(0,i)=sqrt(0.5/correlations(0,i));
  scaler.unScaleXrDist(temp_out_corr_lengths);
  
  std::ostringstream oss;
  oss << "--- Surfpack Kriging Diagnostics ---\n";
  oss << "GEK: #pts="<< numPoints <<"; used " << numEqnKeep << "/" << numEqnAvail << " eqns; Correlation lengths=(" << temp_out_corr_lengths(0,0);
  for(int i=1; i<numVarsr; ++i)
    oss << ", " << temp_out_corr_lengths(0,i);
  oss << "); unadjusted variance=" << estVarianceMLE * scaler.unScaleFactorVarY() << "; log(likelihood)=" << likelihood << "; the trend is a ";
  if(polyOrder>1) {
    if(ifReducedPoly==true)
      oss << "reduced_";
    else oss <<"full ";
  }
  oss << "polynomial of order=" << polyOrder << " (order "<< polyOrderRequested
      << " was desired); rcond(R)=" << rcondR << "; rcond(Gtran_Rinv_G)=" 
      << rcond_Gtran_Rinv_G << "; nugget=" << nug << ".\n";
  
  oss << "Beta= (" << betaHat(0,0);
  //int ntrend=getNTrend(); 
  for(int i=1; i<nTrend; ++i)
    oss << "," << betaHat(i,0);
  oss << ")\n";
  //if(numEqnKeep!=numEqnAvail)
  //for(int i=0; i<numEqnKeep; ++i)
  // oss << "(" << i<< ":" << iptIderKeep(i,0) << "," << iptIderKeep(i,1) << ")\n";
  oss << "------------------------------------\n";    

  return (oss.str());  
}

// BMA TODO: combine these two functions?

/// evaluate (y) the Kriging Model at a single point (xr)
double GradKrigingModel::evaluate(const MtxDbl& xr) const
{
  /*
  double singular_y;
  if(scaler.isYSingular(0,singular_y)) {
    //you wouldn't want to do this for gradient based Kriging
    //if gradients of y were used as inputs
    return singular_y;
  }
  */

  //assert( (numVarsr == xr.getNCols()) && (xr.getNRows() == 1) );
  //int ntrend = getNTrend(); 
  MtxDbl g(1, nTrend), r(1, numRowsR);

  /*
  printf("double evaluate()\n");
  printf("xr=[%20.14g", xr(0,0));
  for(int i=1; i<numVarsr; ++i)
    printf(", %20.14g",xr(0,i));
    printf("]\n");
  */
  

  //MtxDbl xr_scaled(xr);
  if(scaler.isUnScaled()) {
    eval_trend_fn(g, xr);
    correlation_matrix(r, xr);
  }
  else{
    MtxDbl xr_scaled(xr);
    scaler.scaleXrOther(xr_scaled);
    eval_trend_fn(g, xr_scaled);
    correlation_matrix(r, xr_scaled);
  }
      
  if(nug>0.0)
    apply_nugget_eval(r);
  
  //printf("size(xr)={%d %d} size(xr_scaled)={%d %d} size(r)={%d,%d} size(rhs)={%d,%d}\n",
  // xr.getNRows(),xr.getNCols(),
  // xr_scaled.getNRows(),xr_scaled.getNCols(),
  // r.getNRows(),r.getNCols(),
  // rhs.getNRows(),rhs.getNCols());

  double y = dot_product(g, betaHat) + dot_product(r, rhs);

  //double yus=scaler.unScaleYOther(y);
  //printf("] y=%g\n",yus);
  //printf("y=%g yunscaled=%g\n",y,yus);
  //return yus;
  //printf("GEK::evaluate({%g->%g}",xr(0,0),xr_scaled(0,0));
  //for(int ivar=1; ivar<numVarsr; ++ivar)
  //printf(",{%g->%g}",xr(0,ivar),xr_scaled(0,ivar));
  //printf(")=(scaled_y=%g),scaled_y=%g\n",y,scaler.unScaleYOther(y));

  return (scaler.unScaleYOther(y));
}


/// evaluate (y) the Kriging Model at a collection of points (xr)
MtxDbl& GradKrigingModel::evaluate(MtxDbl& y, const MtxDbl& xr) const
{
  int nrowsxr = xr.getNRows();
  //printf("nrowsxr=%d nvarsrxr=%d",nrowsxr,xr.getNCols());

  y.newSize(nrowsxr, 1);
  /*
  double singular_y;
  if(scaler.isYSingular(0,singular_y)) {
    //you wouldn't want to do this for gradient based Kriging
    //if gradients of y were used to build model
    for(int i=0; i<nrowsxr; ++i)
      y(i)=singular_y;
    return y;
  }
  */  

  //assert( (numVarsr == xr.getNCols()) && (xr.getNRows() == 1) );
  //int ntrend = getNTrend(); 
  MtxDbl g(nrowsxr, nTrend), r(nrowsxr, numRowsR);

  MtxDbl xr_scaled(xr); //for debug only, otherwise should be inside "else"
  if(scaler.isUnScaled()) {
    eval_trend_fn(g, xr);
    correlation_matrix(r, xr);
  }
  else{
    //MtxDbl xr_scaled(xr);
    scaler.scaleXrOther(xr_scaled);
    eval_trend_fn(g, xr_scaled);
    correlation_matrix(r, xr_scaled);
  }

  if(nug>0.0)
    apply_nugget_eval(r);
  
  //y=0.0*y+1.0*g*betaHat => y=g*betaHat
  matrix_mult(y, g , betaHat, 0.0, 1.0);
  
  //MtxDbl gbetaHat(y);

  //y=1.0*y+1.0*r*rhs where rhs=R^-1*(Y-G(XR)*betaHat), initial y=g*betaHat => y=g*betaHat+r*rhs
  matrix_mult(y, r, rhs    , 1.0, 1.0);

  /*    
  FILE *fp=fopen("test_eval_scaled_y.txt","w");
  fprintf(fp,"xr_scaled=");
  for(int i=0; i<xr_scaled.getNRows(); ++i) {
    fprintf(fp,"\n%-12.6g",xr_scaled(i,0));
    for(int j=1; j<xr_scaled.getNCols(); ++j)
      fprintf(fp," %-12.6g",xr_scaled(i,j));
  }
  
  fprintf(fp,"\n\ng=");
  for(int i=0; i<g.getNRows(); ++i) {
    fprintf(fp,"\n%-12.6g",g(i,0));
    for(int j=1; j<g.getNCols(); ++j)
      fprintf(fp," %-12.6g",g(i,j));
  }

  fprintf(fp,"\n\nbetaHat=");
  for(int i=0; i<betaHat.getNRows(); ++i)
    fprintf(fp,"\n%-12.6g",betaHat(i,0));

  fprintf(fp,"\n\ngbetaHat=");
  for(int i=0; i<gbetaHat.getNRows(); ++i)
    fprintf(fp,"\n%-12.6g",gbetaHat(i,0));

  fprintf(fp,"\n\nr=");
  for(int i=0; i<r.getNRows(); ++i) {
    fprintf(fp,"\n%-12.6g",r(i,0));
    for(int j=1; j<r.getNCols(); ++j)
      fprintf(fp," %-12.6g",r(i,j));
  }

  fprintf(fp,"\n\nrhs=");
  for(int i=0; i<rhs.getNRows(); ++i)
    fprintf(fp,"\n%-12.6g",rhs(i,0));
  
  fprintf(fp,"\n\nscaled_y=");
  for(int i=0; i<y.getNRows(); ++i)
    fprintf(fp,"\n%-12.6g",y(i,0));
  fflush(fp);
  fclose(fp);
  assert(false);
  */

  scaler.unScaleYOther(y);
  //printf("y is correct for ValidateMain because it isn't being unscaled\n");

  return y;
}

MtxDbl& GradKrigingModel::evaluate_d1y(MtxDbl& d1y, const MtxDbl& xr) const
{
  int nrowsxr = xr.getNRows();
  d1y.newSize(nrowsxr, numVarsr);
  /*
  double singular_y;
  if(scaler.isYSingular(0,singular_y)) {
    //you wouldn't want to do this for gradient based Kriging
    //if gradients of y were used to build model
    d1y.zero();
    return d1y;
  }
  */

  /*
  printf("evaluate_d1y()\n");
  for(int i=0; i<numPoints; ++i) {
    printf("XR(%3d,:)=[%12.6g",i,XR(i,0));
    for(int j=1; j<numVarsr; ++j) 
      printf(", %12.6g",XR(i,j));
    printf("] Y(%3d)=%12.6g\n",i,Y(i,0));
  }
  */

  MtxDbl xr_scaled(xr);  
  if(~(scaler.isUnScaled())) {
    //printf("scaling xr_scaled\n");
    scaler.scaleXrOther(xr_scaled);
  }
  
  /*
  printf("xr       =[%12.6g, %12.6g]\n",xr(0,0),xr(0,1));
  printf("xr_scaled=[%12.6g, %12.6g]\n",xr_scaled(0,0),xr_scaled(0,1));
  */

  //assert(numVarsr == xr.getNCols());
  int nder=num_multi_dim_poly_coef(numVarsr,-1);
  MtxInt der(nder,numVarsr); 
  multi_dim_poly_power(der,numVarsr,-1); //equivalent to der.identity();

  //int ntrend=getNTrend();

  evaluate_poly_der(d1y,Poly,der,betaHat,xr_scaled);
  
#ifdef __KRIGING_DER_TEST__
  assert((d1y.getNRows()==nrowsxr)&&(d1y.getNCols()==numVarsr)&&
	 (der.getNRows()==nder)&&(der.getNCols()==numVarsr));
  MtxDbl d1yalt(nrowsxr*numVarsr,1);
  int npoly=Poly.getNRows();
  MtxDbl d1g(nrowsxr*numVarsr,npoly);
  matrix_mult(d1yalt,evaluate_poly_der_basis(d1g,Poly,der,xr_scaled),betaHat);
  assert((d1yalt.getNElems()==nrowsxr*numVarsr)&&
	 (d1g.getNRows()==nrowsxr*numVarsr)&&(d1g.getNCols()==npoly));
  d1yalt.reshape(nrowsxr,numVarsr);
  for(int ider=0; ider<numVarsr; ++ider)
    for(int ipt=0; ipt<numVarsr; ++ipt) {
      printf("\n d1y(ipt=%d,ider=%d)=%g ",ipt,ider,d1y(ipt,ider));
      if_close_enough(d1y(ipt,ider),d1yalt(ipt,ider));
    }
  printf("\n");
#endif
  
  MtxDbl r(nrowsxr,numRowsR);
  correlation_matrix(r, xr_scaled);
  apply_nugget_eval(r);
  MtxDbl d1r(nrowsxr,numRowsR);
  MtxDbl work(nrowsxr,numPoints); //should this numPoints be NumRowsR?
  MtxDbl temp_vec(nrowsxr,1);


  int ivar;
  for(int ider=0; ider<nder; ++ider) {

    //find the single dimension we are taking the first derviative of
    for(ivar=0; ivar<numVarsr; ++ivar)
      if(der(ider,ivar)>0)
	break;
    //printf("ivar=%d ",ivar);

    double d1y_unscale_factor=scaler.unScaleFactorDerY(ivar);
    //printf("d1y_usf=%g\n",d1y_unscale_factor);

    dcorrelation_matrix_dxI(d1r, r, xr_scaled, work, ider);
    matrix_mult(temp_vec,d1r,rhs);

    for(int ipt=0; ipt<nrowsxr; ++ipt)
      d1y(ipt,ider)=(d1y(ipt,ider)+temp_vec(ipt,0))*d1y_unscale_factor;
  }
  /*
  printf("d1y(0,:)=[%g",d1y(0,0));
  for(int ider=1; ider<numVarsr; ++ider)
    printf(", %g",d1y(0,ider));
  printf("]\n");
  */   
  return d1y;
}

MtxDbl& GradKrigingModel::evaluate_d2y(MtxDbl& d2y, const MtxDbl& xr) const
{
  int nrowsxr=xr.getNRows();
  int nder=num_multi_dim_poly_coef(numVarsr,-2);
  d2y.newSize(nrowsxr,nder);
  /*
  double singular_y;
  if(scaler.isYSingular(0,singular_y)) {
    //you wouldn't want to do this for gradient based Kriging
    //if gradients of y were used as inputs
    d2y.zero();
    return d2y;
  }
  */

  MtxDbl xr_scaled(xr);  
  if(~(scaler.isUnScaled())) 
    scaler.scaleXrOther(xr_scaled);
  //assert(numVarsr == xr.getNCols());

  MtxInt der(nder,numVarsr); 
  MtxInt thisder(1,numVarsr);
  multi_dim_poly_power(der,numVarsr,-2); 

  evaluate_poly_der(d2y,Poly,der,betaHat,xr_scaled);
  
#ifdef __KRIGING_DER_TEST__
  assert((d2y.getNRows()==nrowsxr)&&(d2y.getNCols()==nder)&&
	 (der.getNRows()==nder)&&(der.getNCols()==numVarsr));
  MtxDbl d2yalt(nrowsxr*nder,1);
  int npoly=Poly.getNRows();
  MtxDbl d2g(nrowsxr*nder,npoly);
  matrix_mult(d2yalt,evaluate_poly_der_basis(d2g,Poly,der,xr),betaHat);
  assert((d2yalt.getNElems()==nrowsxr*nder)&&
	 (d2g.getNRows()==nrowsxr*nder)&&(d2g.getNCols()==npoly));
  d2yalt.reshape(nrowsxr,nder);
  for(int ider=0; ider<nder; ++ider)
    for(int ipt=0; ipt<numVarsr; ++ipt) {
      printf("\n d2y(ipt=%d,ider=%d)=%g ",ipt,ider,d2y(ipt,ider));
      int ifok=if_close_enough(d2y(ipt,ider),d2yalt(ipt,ider));
    }
  printf("\n");
#endif

  MtxDbl r(nrowsxr,numRowsR);
  correlation_matrix(r, xr);
  apply_nugget_eval(r);
  MtxDbl d1r(nrowsxr,numRowsR);
  MtxDbl d2r(nrowsxr,numRowsR);
  MtxDbl temp_vec(nrowsxr,1);
  MtxDbl work(nrowsxr,numPoints);

  int ivar, ivarold=-1, jvar;
  for(int ider=0; ider<nder; ++ider) {

    der.getRows(thisder,ider);
    double d2y_unscale_factor=scaler.unScaleFactorDerY(thisder);

    //find the first dimension we are taking a first derviative of
    for(ivar=0; ivar<numVarsr; ++ivar)
      if(der(ider,ivar)>0)
	break;

    if(ivar!=ivarold) {
      ivarold=ivar;
      dcorrelation_matrix_dxI(d1r, r, xr_scaled, work, ivar);
    }

    //find the second dimension we are taking a first derivative of
    if(der(ider,ivar)==2)
      jvar=ivar;
    else
      for(jvar=ivar+1; jvar<numVarsr; ++jvar)
	if(der(ider,jvar)>0)
	  break;
    
    d2correlation_matrix_dxIdxK(d2r, d1r, r, xr_scaled, work, ivar, jvar);
    
    matrix_mult(temp_vec,d2r,rhs);

    for(int ipt=0; ipt<nrowsxr; ++ipt)
      d2y(ipt,ider)=(d2y(ipt,ider)+temp_vec(ipt,0))*d2y_unscale_factor;
  }

  return d2y;
}



/// matrix Ops evaluation of adjusted variance at a single point
double GradKrigingModel::eval_variance(const MtxDbl& xr) const
{
  //MtxDbl m_adj_var(1,1);
  //eval_variance(m_adj_var, xr);
  //double adj_var=m_adj_var(0,0);

  /*
  double singular_y;
  if(scaler.isYSingular(0,singular_y)) {
    //you wouldn't want to do this for gradient based Kriging
    //if gradients of y were used as inputs
    //printf("NKM eval_variance: y is singular\n");
    return 0;
  }
  */

  

  //assert( (numVarsr == xr.getNCols()) && (xr.getNRows() == 1) );
  //int ntrend = getNTrend(); 
  MtxDbl g_minus_r_Rinv_G(1, nTrend), r(1, numRowsR);

  MtxDbl xr_scaled(xr);
  if(scaler.isUnScaled()) {
    //printf("unscaled path\n");
    eval_trend_fn(g_minus_r_Rinv_G, xr);
    correlation_matrix(r, xr);
  }
  else{
    //prinf("scaled path\n");
    //MtxDbl xr_scaled(xr);
    scaler.scaleXrOther(xr_scaled);
    eval_trend_fn(g_minus_r_Rinv_G, xr_scaled);
    correlation_matrix(r, xr_scaled);
  }

  MtxDbl tempa(numRowsR,1);
  MtxDbl tempb(nTrend,1);

  if(nug>0.0) 
    apply_nugget_eval(r);

  //matrix_mult(tempa,Rinv,r,0.0,1.0,'N','T');
  solve_after_Chol_fact(tempa,RChol,r,'T');
  
  matrix_mult(g_minus_r_Rinv_G,r,Rinv_G,1.0,-1.0);
  //matrix_mult(tempb,Gtran_Rinv_G_inv,g_minus_r_Rinv_G,0.0,1.0,'N','T');
  solve_after_Chol_fact(tempb,Gtran_Rinv_G_Chol,g_minus_r_Rinv_G,'T');

  

  double unscale_factor_vary=scaler.unScaleFactorVarY();

  double unadj_var=estVarianceMLE*unscale_factor_vary;
  
  double dtemp1=dot_product(tempa,r);
  double dtemp2=dot_product(tempb,g_minus_r_Rinv_G);
  double adj_var=unadj_var*(1.0-dtemp1+dtemp2);
  //if(!(adj_var>0.0)) {
  //printf("adj_var=%g estVarianceMLE=%g rcondR=%g unscale_factor_vary=%g\n",adj_var,estVarianceMLE,rcondR,unscale_factor_vary); 
  //fflush(stdout);
  //}
  
  adj_var=std::fabs(adj_var);
  if(adj_var<0.0) {
    /*
    printf("NKM GEK setting adj_var to zero adj_var=%g unadj_var=%g rcondR=%g\n",adj_var,unadj_var,rcondR); 
    printf("isUnScaled=%d\ncorr_len=[%12.6g",(scaler.isUnScaled())?1:0,1.0/std::sqrt(2.0*correlations(0,0)));
    for(int i=1; i<numVarsr; ++i)
      printf(" %12.6g",1.0/std::sqrt(2.0*correlations(0,i)));
    printf("]\nnrowsxr=%d\n",xr_scaled.getNRows());
    printf("xr=|%12.6g",xr_scaled(0,0));
    for(int j=1; j<numVarsr; ++j)
      printf(" %12.6g",xr_scaled(0,j));
    printf("|\nnumPoints=%d numEqnKeep=%d\n(ipt)|XR|\n",numPoints,numEqnKeep);
    for(int i=0; i<numPoints; ++i) {
      int ipt=iEqnKeep(i*(1+numVarsr),0);
      printf("(%4d)|%12.6g",ipt,XR(ipt,0));
      for(int j=1; j<numVarsr; ++j)
	printf(" %12.6g",XR(ipt,j));
      printf("|\n");
    }
    printf("rU=[...\n");
    for(int i=0; i<numEqnKeep; ++i){
      printf("%12.6g *****",r(i,0));
      for(int j=0; j<i; ++j)
	printf(" %12.6g",0.0);
      for(int j=i; j<numEqnKeep; ++j)
	printf(" %12.6g",RChol(i,j));
      printf(";...\n");
    }
    printf("];\n");

    printf("unscale_factor_vary=%g nTrend=%d numRowsR=%d dtemp1=%g dtemp2=%g\n",unscale_factor_vary,nTrend,numRowsR,dtemp1,dtemp2);
    assert(false);
    */
    adj_var=0.0;
  }
  else if(adj_var==0.0)
    printf("NKM GEK adj_var is zero =%g\n",adj_var);
  else if(!(adj_var>=0.0))
    printf("double NKM_GradKrigingModel::eval_variance(...) adj_var=nan rcondR=%g\n",rcondR);

  return adj_var;
}

/// matrix Ops (as much as possible with BLAS and LAPACK) evaluation of adjusted variance for a collection of points... The MATLAB would be estVarianceMLE*(1-sum((r/R).*r,2)+sum((g_minus_r_Rinv_G/(Gtran_Rinv_G)).*g_minus_r_Rinv_G,2) unfortunately there's not a convenient way to do it with BLAS & LAPACK
MtxDbl& GradKrigingModel:: eval_variance(MtxDbl& adj_var, const MtxDbl& xr) const
{
  int nrowsxr=xr.getNRows();
  adj_var.newSize(nrowsxr,1);

  /*
  double singular_y;
  if(scaler.isYSingular(0,singular_y)) {
    //you wouldn't want to do this for gradient based Kriging
    //if gradients of y were used as inputs
    adj_var.zero();
    return adj_var;
  }
  */

  //assert( (numVarsr == xr.getNCols()) && (xr.getNRows() == 1) );
  //int ntrend = getNTrend(); 
  MtxDbl g_minus_r_Rinv_G(nrowsxr, nTrend), r(nrowsxr, numRowsR);

  if(scaler.isUnScaled()) {
    eval_trend_fn(g_minus_r_Rinv_G, xr);
    correlation_matrix(r, xr);
  }
  else{
    MtxDbl xr_scaled(xr);
    scaler.scaleXrOther(xr_scaled);
    eval_trend_fn(g_minus_r_Rinv_G, xr_scaled);
    correlation_matrix(r, xr_scaled);
  }

  int i,j;
  MtxDbl tempa(numRowsR,nrowsxr); 
  MtxDbl tempb(nTrend,nrowsxr);
  double var_unscale_factor=scaler.unScaleFactorVarY();

  if(nug>0.0)
    apply_nugget_eval(r);

  //matrix_mult(tempa,Rinv,r,0.0,1.0,'N','T');  
  solve_after_Chol_fact(tempa,RChol,r,'T');

  matrix_mult(g_minus_r_Rinv_G,r,Rinv_G,1.0,-1.0);
  //matrix_mult(tempb,Gtran_Rinv_G_inv,g_minus_r_Rinv_G,0.0,1.0,'N','T');
  solve_after_Chol_fact(tempb,Gtran_Rinv_G_Chol,g_minus_r_Rinv_G,'T');

  for(i=0; i<nrowsxr; ++i) {
    //saved 2*nrowsxr loops
    adj_var(i,0)=1.0-r(i,0)*tempa(0,i)+g_minus_r_Rinv_G(i,0)*tempb(0,i);

    for(j=1; j<numRowsR; ++j)
      adj_var(i,0)-=r(i,j)*tempa(j,i); //looks a lot like matrix mult but only N^2 ops

    for(j=1; j<nTrend; ++j)
      adj_var(i,0)+=g_minus_r_Rinv_G(i,j)*tempb(j,i); //looks a lot like matrix mult but only N^2 ops

    adj_var(i,0)*=estVarianceMLE*var_unscale_factor;
  }

  for(i=0; i<nrowsxr; ++i)
    if(adj_var(i,0)<0.0)
      adj_var(i,0)=0.0;
    else if(!(adj_var(i,0)>=0.0))
      printf("MtxDbl& NKM_GradKrigingModel::eval_variance(...) adj_var(%d)=nan rcondR=%g\n",i,rcondR);
	

  return adj_var;
}


/*
VecDbl GradKrigingModel::gradient(const VecDbl& x) const
{
  assert(!x.empty());
  assert(x.size()+1==betaHat.size()); //true for linear trend function; KRD added
  std::cout << "IN gradient x[0] = " << x[0] << std::endl;
  assert(rhs.size() == bs.centers.size());
  VecUns diff_var(1,0); // variable with which to differentiate
  VecDbl result(x.size(),0.0);
  for (unsigned i = 0; i < x.size(); i++) {
    diff_var[0] = i;
    result[i] = betaHat[i+1]; //true for linear trend function; KRD added
    for (unsigned j = 0; j < bs.centers.size(); j++) {
      result[i] += rhs[j]*bs.deriv(j,x,diff_var);
    }
  }
  return result;
}

std::string GradKrigingModel::asString() const
{
  std::ostringstream os;
  os << "\ncenters:\n" << bs.asString() << "rhs: ";
  copy(rhs.begin(),rhs.end(),std::ostream_iterator<double>(os," "));
  os << "\n";
  return os.str();
}
*/


/** this function scales the input matrix, r, in place by 1.0/(1.0+nug), 
    where nug is the nugget (modifying the correlation matrix by the 
    inclusion of a nugget causes the GradKrigingModel to smooth data, i.e. 
    approximate, rather than interpolate, it can also be used to fix a 
    poorly conditioned correlation matrix).  Typically the input matrix 
    will be the correlation matrix, r, for arbitrary points at which to 
    EVALuate the model (hence the name), but scaling the _derivative_ of R 
    (the correlation matrix for data the model was built from) with respect 
    to theta is also common, since the additive term, nug*I, is constant with 
    respect to theta (i.e. the additive term for the build data correlation 
    matrix R drops out of its derivative).  The convention is that capital
    matrices are for the data the model is built from, lower case matrices 
    are for arbitrary points to evaluate the model at. */
MtxDbl& GradKrigingModel::apply_nugget_eval(MtxDbl& r) const {
  if(!(nug>0.0))
    return r;

  //printf("apply_nugget_eval\n");

  //int nelem=r.getNElems();
  int nrowsr=r.getNRows();
  int ncolsr=r.getNCols();
  double temp_dbl=1.0/(1.0+nug);

  //for(int ij=0; ij<nelem; ++ij) 
  //r(ij)*=temp_dbl;
  for(int j=0; j<ncolsr; ++j)
    for(int i=0; i<nrowsr; ++i)
      r(i,j)*=temp_dbl;

  //printf("apply_nugget_eval temp=%g\n",temp);

  return r;
}

/** set R=(1.0+nug)^-1*(R+nug*I), where R is the correlation matrix for the
    data that the model is built from.  This function preserved the original 
    R.  Modifying the correlation matrix by the inclusion of a nugget causes 
    the GradKrigingModel to smooth the data, i.e. approximate it rather than 
    interpolate it, it can also be used to fix an ill conditioned correlation 
    matrix.  The convention is that capital matrices are for data the model 
    is built from, lower case matrices are for arbitrary points to evaluate 
    the model at */
void GradKrigingModel::apply_nugget_build() {
  if(!(nug>0.0)) return;
  
  //int nrowsR=R.getNRows();
  //assert(nrowsR==R.getNCols());
  //int nelemsR=numRowsR*numRowsR;

  double temp_dbl=1.0/(1.0+nug);
  /*
  int ij;
  for(ij=0; ij<nelemsR; ++ij)
    R(ij)*=temp_dbl;
  */
  for(int j=0; j<numRowsR; ++j)
    for(int i=0; i<numRowsR; ++i)
      R(i,j)*=temp_dbl;

  //the "paranoid" part of my mind wonders if there would be less round off
  //error if I added the nugget to the diagonal BEFORE scaling, the pragmatic
  //part of my mind says it shouldn't matter and doing it this way is faster
  temp_dbl*=nug;
  for(int j=0; j<numPoints; ++j) //numPoints is correct because derivatives 
    //(i.e. j >= numPoints) of a constant is zero so don't add it there
    R(j,j)+=temp_dbl; 
    
  return;
}



/** r (lower case r) is the correlation matrix between the
    interpolation points and data points, it used to EVALUATE but not
    construct the emulator's Gaussian process error model
    i.e. E(y(xr)|Y(XR))=g(xr)*betaHat+r*R^-1*eps where
    eps=(Y-G(XR)*betaHat), to be more specific
    r(i,j)=r(xr(i),XR(j))=exp(sum_k -theta(k)*(xr(i,k)-XR(j,k))^2) KRD
    wrote this */
MtxDbl& GradKrigingModel::correlation_matrix(MtxDbl& r, const MtxDbl& xr) const
{
  //int nrowsXR=XR.getNRows(); //data points used to build model
  int nrowsxr=xr.getNRows(); //points at which we are evalutating the model
  //assert(xr.getNCols()==numVarsr);
  r.newSize(nrowsxr,numEqnKeep);

  int i,j,k, Jder;

  double temp_double;
  
  int j_last_func_val=0;
  int jpt;
  if(numVarsr==1) {
    for(j=0; j<numEqnKeep; ++j) {
      jpt=iptIderKeep(j,0);
      Jder=iptIderKeep(j,1);
      //printf("corr_mtx(r,xr)::jpt=%d Jder=%d\n",jpt,Jder);
      
      if(Jder==-1) {
	//printf("corr_mtx(r,xr):: case 1\n");
	//correlation between function values
	j_last_func_val=j;
	for(i=0; i<nrowsxr; ++i) {
	  temp_double=xr(i,0)-XR(jpt,0);
	  r(i,j)=std::exp(-correlations(0,0)*temp_double*temp_double); 
	}
      }
      else{
	//printf("corr_mtx(r,xr):: case 2\n");
	//correlation between the prediction point's function values and the
	//"Jder" derivative of the build point
	//relies on knowing that if a derivative equation occurs, the previous
	//function value will be the one it's a derivative of
	temp_double=2.0*correlations(0,Jder);
	for(i=0; i<nrowsxr; ++i)
	  r(i,j)=temp_double*(xr(i,Jder)-XR(jpt,Jder))*r(i,j_last_func_val);
      }
    }
  }
  else{
    for(j=0; j<numEqnKeep; ++j) {
      jpt=iptIderKeep(j,0);
      Jder=iptIderKeep(j,1);
      //printf("corr_mtx(r,xr)::jpt=%d Jder=%d\n",jpt,Jder);
      
      if(Jder==-1) {
	//printf("corr_mtx(r,xr):: case 1\n");
	//correlation between function values
	j_last_func_val=j;
	for(i=0; i<nrowsxr; ++i) {
	  temp_double=xr(i,0)-XR(jpt,0);
	  r(i,j)=-correlations(0,0)*temp_double*temp_double; //=- is correct
	}
	for(k=1; k<numVarsr-1; ++k)
	  for(i=0; i<nrowsxr; ++i) {
	    temp_double=xr(i,k)-XR(jpt,k);
	    r( i, j)-=correlations(0,k)*temp_double*temp_double; //-= is correct
	  }
	k=numVarsr-1;
	for(i=0; i<nrowsxr; ++i) {
	  temp_double=xr(i,k)-XR(jpt,k);
	  r( i, j)=std::exp(r( i, j)-correlations(0,k)*temp_double*temp_double);
	}
      }
      else{
	//printf("corr_mtx(r,xr):: case 2\n");
	//correlation between the prediction point's function values and the
	//"Jder" derivative of the build point
	//relies on knowing that if a derivative equation occurs, the previous
	//function value will be the one it's a derivative of
	temp_double=2.0*correlations(0,Jder);
	for(i=0; i<nrowsxr; ++i)
	  r(i,j)=temp_double*(xr(i,Jder)-XR(jpt,Jder))*r(i,j_last_func_val);
      }
    }
  }

  return r;
}


MtxDbl& GradKrigingModel::dcorrelation_matrix_dxI(MtxDbl& dr, const MtxDbl& r, 
						  const MtxDbl& xr, 
						  MtxDbl& workI, int Ider) const
{
  int nrowsXR=XR.getNRows(); //data points used to build model
  int nrowsxr=xr.getNRows(); //points at which we are evalutating the model
  //assert(xr.getNCols()==numVarsr);
  dr.newSize(nrowsxr,numEqnKeep);
  workI.newSize(nrowsxr,nrowsXR);

  double twoThetaI=2.0*correlations(0,Ider);
  //now I don't need to touch correlations for the rest of the function call

  for(int jpt=0; jpt<nrowsXR; ++jpt)
    for(int i=0; i<nrowsxr; ++i)
      workI(i,jpt)=-twoThetaI*(xr(i,Ider)-XR(jpt,Ider));
  //now I don't need to touch xr or XR for the rest of the function call 
  //(because I have workI) it's possible that adding workI slowed things down,
  //we now need only 1 memory access instead of 2 but need nrowsxr*nrowsXR 
  //memory instead of nrowsxr+nrowsXR the good news that the nrowsxr*nrowsXR 
  //is laid out linearly in memory in the order it will be needed so the
  //compiler should be able to optimize that pretty well

  int j_last_func_val;
  for(int j=0; j<numEqnKeep; ++j) {
    int jpt=iptIderKeep(j,0);
    int Jder=iptIderKeep(j,1);

    if(Jder==-1) {
      j_last_func_val=j;
      //first derivative with respect to Ider component (variable) of first 
      //input, xr(:,Ider), no derivatives W.R.T. the second input (XR)
      for(int i=0; i<nrowsxr; ++i)
	dr(i,j)=workI(i,jpt)*r(i,j);
    }
    else{
      //now do the mixed second derivatives submatrices  
      if(Jder==Ider) {
	//relies on knowing that if a derivative W.R.T. XR equation occurs, 
	//the previous function value will be the one it's a derivative of
	for(int i=0; i<nrowsxr; ++i) 
	  dr(i,j)=workI(i,jpt)*r(i,j)+twoThetaI*r(i,j_last_func_val);
      }
      else
	for(int i=0; i<nrowsxr; ++i) 
	  dr(i,j)=workI(i,jpt)*r(i,j);
    }
  }
  
  return dr;
}

MtxDbl& GradKrigingModel::d2correlation_matrix_dxIdxK(MtxDbl& d2r, const MtxDbl& drI, const MtxDbl& r, const MtxDbl& xr, MtxDbl& workK, int Ider, int Kder) const
{
  int nrowsXR=XR.getNRows(); //data points used to build model
  int nrowsxr=xr.getNRows(); //points at which we are evalutating the model
  d2r.newSize(nrowsxr,numEqnKeep);
  workK.newSize(nrowsxr,nrowsXR);

  double twoThetaK=2.0*correlations(0,Kder);
  //now I don't need to touch correlations for the remainder of the function 
  //call

  for(int jpt=0; jpt<nrowsXR; ++jpt)
    for(int i=0; i<nrowsxr; ++i)
      workK(i,jpt)=-twoThetaK*(xr(i,Kder)-XR(jpt,Kder));
  //now I don't need to touch xr or XR for the remainder of the function call


  int j_last_func_val;
  if(Ider==Kder) {
    //second derivative with respect to Ider component (variable) of first 
    //input, d^2[...]/dxr(:,Ider)dxr(:,Kder)  and Ider==Kder
    for(int j=0; j<numEqnKeep; ++j) {
      int jpt=iptIderKeep(j,0);
      int Jder=iptIderKeep(j,1);

      if(Jder==-1) {
	j_last_func_val=j;
	//no derivatives with respect to the second input (XR)
	for(int i=0; i<nrowsxr; ++i)
	  d2r(i,j)=workK(i,jpt)*drI(i,j)-twoThetaK*r(i,j);
      }
      else{
	//now do the mixed third derivatives submatrices (with a first 
	//derivative with respect to XR) so 
	//d^3[r(xr(i,:),XR(jpt,:))]/dXR(jpt,Jder)dxr(i,Ider)dxr(i,Kder) 
	//& Ider==Kder
	if(Jder==Kder) {
	  //relies on knowing that if a derivative W.R.T. XR equation occurs, 
	  //the previous function value will be the one it's a derivative of
	  for(int i=0; i<nrowsxr; ++i)
	    d2r(i,j)=workK(i,jpt)*drI(i,j)-twoThetaK*r(i,j)
		+twoThetaK*drI(i,j_last_func_val);
	}
	else
	  for(int i=0; i<nrowsxr; ++i)
	    d2r(i,j)=workK(i,jpt)*drI(i,j)-twoThetaK*r(i,j);
      }
    }
  }
  else{
    //mixed second derivatives with respect to Ider and Kder components 
    //(variables) of the first input,d^2[...]/dxr(:,Ider)dxr(:,Kder) 
    //Ider =/= Kder
    for(int j=0; j<numEqnKeep; ++j) {
      int jpt=iptIderKeep(j,0);
      int Jder=iptIderKeep(j,1);
      
      if(Jder==-1) {
	j_last_func_val=j;
	//no derivatives with respect to the second input (XR)
	for(int i=0; i<nrowsxr; ++i)
	  d2r(i,j)=workK(i,jpt)*drI(i,j);
      }
      else{
	//now do the mixed third derivatives submatrices (with a first 
	//derivative with respect to XR) so 
	//d^3[r(xr(i,:),XR(jpt,:))]/dXR(jpt,Jder)dxr(i,Ider)dxr(i,Kder) 
	//& Ider=/=Kder
	if(Jder==Kder) {
	  //relies on knowing that if a derivative W.R.T. XR equation occurs, 
	  //the previous function value will be the one it's a derivative of
	  for(int i=0; i<nrowsxr; ++i)
	    d2r(i,j)=workK(i,jpt)*drI(i,j)
	      +twoThetaK*drI(i,j_last_func_val);
	}
	else
	  for(int i=0; i<nrowsxr; ++i)
	    d2r(i,j)=workK(i,jpt)*drI(i,j);		
      }
    }
  }
  
  return d2r;
}

/** this function is typically used during emulator construction, the below
    the diagonal portion of R = exp(Z*theta), where R is symmetric with 1's 
    on the diagonal, theta is the vector of correlations and the Z matrix is 
    defined as Z(ij,k)=-(XR(i,k)-XR(j,k))^2 where ij counts down columns 
    from the element below the diagonal and continues from one column to the 
    next, Z*theta is matrix vector multiplication to be performed efficiently 
    by BLAS, V=Z*theta is a vector with nchoosek(numPoints,2) elements.  We 
    need to copy exp(V(ij)) to R(i,j) and R(j,i) to produce R. The Z matrix 
    is produced by GradKrigingModel::gen_Z_matrix()
 */
void GradKrigingModel::correlation_matrix(const MtxDbl& theta)
{
  //correlations.copy(theta);
  int nrowsZ=Z.getNRows();
  //printf("nrowsZ=%d numRowsR=%d\n",nrowsZ,numRowsR);

  Ztheta.newSize(nrowsZ,1);
  matrix_mult(Ztheta,Z,theta,0.0,1.0,'N','T');

  R.newSize(numEqnAvail,numEqnAvail); 
  //the convention is that each index into R has 2 characters R(__,__) 
  //eg R(Ii,Jj). The second charater is lower case and indicates the "point"
  //index INTO a derivative submatrix.  The first character in each index is 
  //upper case and indicates the derivative submatrix itself, more 
  //specifically which input variable the derivative is with respect to. 
  //For example R(Ii,Jj) indicates the ith row (point) and jth column (point) 
  //of the (I_,J_) submatrix which means 
  //d^2[r(xr(i,:),XR(j,:))]/dXR(i,I)dXR(j,J), 
  //A blank space for the first character of an index means NO 
  //derivative was taken with respect to that input, for example 
  //R( i, j)=r(xr(i,:),XR(j,:))  
  //R(Ii, j)=dr(xr(i,:),XR(j,:))/dXR(i,I)
  //R(Ji,Jj) indicates this is a diagonal submatrix (i.e. it's a mixed 
  //second derivative with respect to the same component of the first 
  //and second input variables)
  //R(Ji,Jj)=d^2[r(xr(i,:),XR(j,:))]/dXR(i,J)dXR(j,J) 
  //which has the extra /additive term, +2*theta(J)*r( i, j) 
  //relative to other mixed second derivatives

  double temp_double; //a temporary variable from which the matrix is filled
  //so I don't occur a read-access performance loss due to copying from one 
  //ELEMENT of the matrix to another, the compiler should know enough to 
  //make "temp_double" a register variable
  int zij=0; //index into Z and twoXRminusXRtran (only the lower 
  //triangular part of these matrices were stored to save space)
  int j; //j "point" index into a submatrix
  //int i; //i "point" index into a submatrix

  //printf("correlation_matrix(): numEqnAvail=%d numPoints=%d\n",numEqnAvail,numPoints);
  //the zeroth order derivative (function value itself) submatrix equals 
  //the correlation matrix without derivatives and is symmetric
  for(j=0; j<numPoints-1; ++j) {
    R( j, j)=1.0;  
    for(int i=j+1; i<numPoints; ++i) {
      temp_double=std::exp(Ztheta(zij,0));
      R( i, j)=temp_double;
      R( j, i)=temp_double;
      //printf("i=%d j=%d zij=%d Ztheta=%g R=%g\n",i,j,zij,Ztheta(zij),R(i,j));
      ++zij;
    }
  }
  j=numPoints-1;
  R( j, j)=1.0;

  //now handle the first order derivative submatrices, indiviually the first
  //order derivative SUBmatrices are anti-symmetric but the whole matrix is
  //symmetric
  int Ii, Ij, Jj, Ji; //first letter identifies index OF derivative submatrix 
  //second letter identifies index INTO derivative SUBmatrix
  for(int Ider=0; Ider<numVarsr; ++Ider) {
    zij=0;
    for(j=0; j<numPoints-1; ++j) {//j<numPoints-1 avoids an i loop of length 0
      //diagonal (_j,_j) of off diagonal (I_, _) submatrix 
      Ij=(Ider+1)*numPoints+j;
      R(Ij, j)=0.0;
      R( j,Ij)=0.0;
      //Ij=(Ider+1)*numPoints+j;
      for(int i=j+1; i<numPoints; ++i) {
	//off diagonal (_i,_j) of off-diagonal (I_, _) submatrix 
	Ii=(Ider+1)*numPoints+i;
	temp_double=-theta(0,Ider)*twoXRminusXRtran(zij,Ider)*R( i, j);
	R(Ii, j)= temp_double; 
	R( j,Ii)= temp_double; //whole R matrix is symmetric
	R(Ij, i)=-temp_double; 
	R( i,Ij)=-temp_double; //off-diagonal 1st order (actually all odd 
	//order) derivative SUBmatrices are anti-symmetric
	++zij;
      }
    }
    //diagonal (_j,_j) of off diagonal (I_, _) submatrix 
    j=numPoints-1; //avoids an i loop of length 0
    Ij=(Ider+1)*numPoints+j;
    R(Ij,j)=0.0;
    R(j,Ij)=0.0;    
  }

  //note that all 2nd order (actually all even order) derivative SUBmatrices 
  //are symmetric because the hadamard product of 2 (actually any even 
  //number of) anti-symmetric matrices is a symmetric matrix
  double two_theta_Jder;
  for(int Jder=0; Jder<numVarsr; ++Jder) {
    //do the on diagonal (J_,J_) submatrix 
    two_theta_Jder=2.0*theta(0,Jder);
    zij=0;
    for(j=0; j<numPoints-1; ++j) { //j<numPoints-1 avoids an i loop of length 0
      //diagonal (_j,_j) of on diagonal (J_,J_) submatrix
      Jj=(Jder+1)*numPoints+j; 
      R(Jj,Jj)=two_theta_Jder; //R(Jj,Jj)=2*theta(0,Jder)*R(j,j); R(j,j)=1; 
      for(int i=j+1; i<numPoints; ++i) {
	//off diagonal (_i,_j) of on-diagonal (J_,J_) submatrix
	Ji=(Jder+1)*numPoints+i;
	temp_double=theta(0,Jder)*twoXRminusXRtran(zij,Jder)*R(Ji, j)+
	  two_theta_Jder*R( i, j);
	R(Ji,Jj)=temp_double;
	R(Jj,Ji)=temp_double;
	++zij;
      }
    }
    //diagonal (_j,_j) of on diagonal (J_,J_) submatrix 
    j=numPoints-1; //avoids an i loop of length 0
    Jj=(Jder+1)*numPoints+j;
    R(Jj,Jj)=two_theta_Jder; //R(j,j)=1 R(Jj,Jj)=2*theta(0,Jder)*R(j,j)


    //do the off diagonal (I_,J_) submatrices
    for(int Ider=Jder+1; Ider<numVarsr; ++Ider) {
      //off diagonal (I_,J_) submatrix
      zij=0;
      for(j=0; j<numPoints-1; ++j) {//j<numPoints-1 avoids an i loop of length 0
	//diagonal (_j,_j) of off-diagonal (I_,J_) submatrix
	Jj=(Jder+1)*numPoints+j; 
	Ij=(Ider+1)*numPoints+j;
	R(Ij,Jj)=0.0;
	R(Jj,Ij)=0.0;


	for(int i=j+1; i<numPoints; ++i) {
	  //off diagonal (_i,_j) of off-diagonal (I_,J_) submatrix 
	  Ii=(Ider+1)*numPoints+i;
	  Ji=(Jder+1)*numPoints+i;
	  temp_double=theta(0,Jder)*twoXRminusXRtran(zij,Jder)*R(Ii, j);
	  R(Ii,Jj)= temp_double; 
	  R(Ij,Ji)= temp_double; 
	  R(Ji,Ij)= temp_double; 
	  R(Jj,Ii)= temp_double; 
	  ++zij;
	}
      }
      //diagonal (_j,_j) of off-diagonal (I_,J_) submatrix
      j=numPoints-1; //avoids an i loop of length 0
      Ij=(Ider+1)*numPoints+j;
      Jj=(Jder+1)*numPoints+j;
      R(Ij,Jj)=0.0;
      R(Jj,Ij)=0.0;
    }
  }

  /*  
  FILE *fp=fopen("gkm_Rmat_check.txt","w");
  fprintf(fp,"theta={%g",theta(0,0));
  for(int i=1; i<numVarsr; ++i)
    fprintf(fp," %g",theta(0,i));
  fprintf(fp,"}\nR=\n");
  for(int i=0; i<numPoints; ++i) {
    fprintf(fp,"%-12.6g", R(i,0));
    //for(int j=1; j<numRowsR; ++j) 
    for(int j=1; j<numEqnAvail; ++j) 
      fprintf(fp," %-12.6g", R(i,j));     
    fprintf(fp,"\n");
  }
  fclose(fp);
  */


#ifdef _FOR_DEBUG_DEVEL_ONLY_
  //for debug only
  MtxDbl rr(numRowsR,numRowsR); //sizing is bad now
  MtxDbl r(numPoints,numRowsR);
  MtxDbl dr(numPoints,numRowsR);
  MtxInt irows(numPoints);
  MtxDbl work(numPoints,numPoints);

  correlations.copy(theta);
  correlation_matrix(r, XR); //this would be a problematic comparison since iptIderKeep should only be defined for the "anchor points", currently that's defaulting to an empty matrix (i still need to add a copy in for the iAnchorPoints array to the constructor)

  for(int i=0; i<numPoints; ++i)
    irows(i,0)=i;
  rr.putRows(r,irows);

  for(int Ider=0; Ider<numVarsr; ++Ider) {
    irows(0,0)=(Ider+1)*numPoints;
    for(int i=1; i<numPoints; ++i)
      irows(i,0)=irows(i-1,0)+1;
    dcorrelation_matrix_dxI(dr,r,XR,work,Ider);
    rr.putRows(dr,irows);
  }

  
  FILE* fp=fopen("CheckGKE_Rmat.txt","w");
  fprintf(fp,"%d\n%d\n\n",numPoints,numVarsr);
  for(int i; i<numVarsr; ++i)
    fprintf(fp,"%.16g ",theta(0,i));
  fprintf(fp,"\n");
  

  for(int i=0; i<numPoints; ++i) {
    fprintf(fp,"\n%.16g",XR(i,0));
    for(j=1; j<numVarsr; ++j)
      fprintf(fp," %.16g",XR(i,j));
  }
  fprintf(fp,"\n");

  for(j=0; j<numVarsr; ++j) {
    fprintf(fp,"\n%.16g",Z(0,j));
    for(int zij=1; zij<nrowsZ; ++zij) 
      fprintf(fp,"\n%.16g",Z(zij,j));
  }
  fprintf(fp,"\n");

  for(j=0; j<numVarsr; ++j) {
    fprintf(fp,"\n%.16g",twoXRminusXRtran(0,j));
    for(int zij=1; zij<nrowsZ; ++zij) 
      fprintf(fp,"\n%.16g",twoXRminusXRtran(zij,j));
  }
  fprintf(fp,"\n");

  /*
  for(int i=0; i<numRowsR; ++i) {
    fprintf(fp,"\n%.16g",R(i,0)); 
    for(j=1; j<numRowsR; ++j) {
          fprintf(fp," %.16g",R(i,j)); 
    }
  }
  fprintf(fp,"\n");
  */

  for(int i=0; i<numEqnAvail; ++i) {
    fprintf(fp,"\n%.16g",R(i,0)); 
    for(j=1; j<numEqnAvail; ++j) {
          fprintf(fp," %.16g",R(i,j)); 
    }
  }
  fprintf(fp,"\n");

  for(int i=0; i<numRowsR; ++i) {
    fprintf(fp,"\n%.16g",rr(i,0)); 
    for(j=1; j<numRowsR; ++j) {
          fprintf(fp," %.16g",rr(i,j)); 
    }
  }
  fflush(fp);
  fclose(fp);
  assert(false);
#endif  
  return; 
}

/** the Z matrix is defined as Z(ij,k)=-(XR(i,k)-XR(j,k))^2 where
    ij=i+j*XR.getNRows(), it enables the efficient repeated calculation
    of the R matrix during model construction:
    R=reshape(exp(Z*theta),XR.getNRows(),XR.getNRows()) where theta is
    the vector of correlations and * is matrix vector multiplication,
    note that the Z matrix is independent of the correlation vector so
    it can be formed once and later during the search for a good
    correlation vector, the matrix vector product Z*theta can be
    performed efficiently (for each member of the set of candidate
    theta vectors) by calling BLAS. Z and XR are member variables so  
    they don't need to be passed in, KRD wrote this,  */
MtxDbl& GradKrigingModel::gen_Z_matrix()
{
  int nrowsXR=XR.getNRows();
  int ncolsXR=XR.getNCols();
  int nrowsZ=nchoosek(nrowsXR,2);

  twoXRminusXRtran.newSize(nrowsZ,ncolsXR);
  Z.newSize(nrowsZ,ncolsXR);


  register double mult_term;
  register int ijk=0;
  double *twoXRminusXRtran_ptr=twoXRminusXRtran.ptr(0,0);
  double *Z_ptr=Z.ptr(0,0); //done for speed
  const double *XR_k_ptr; //done for speed
  for(int k=0; k<ncolsXR; k++) {
    XR_k_ptr=XR.ptr(0,k);
    for(int j=0; j<nrowsXR-1; j++)
      for(int i=j+1; i<nrowsXR; i++) {
	mult_term=XR_k_ptr[i]-XR_k_ptr[j];
	twoXRminusXRtran_ptr[ijk]=2*mult_term;
	Z_ptr[ijk++]=-mult_term*mult_term;
      }
  }
  
  return Z;
}

// BMA TODO: combine shared code from these various functions?
// (consider initializing class members R1 and R2 or something once
// per eval and reusing)

// BMA TODO: add code for likelihood function calculation

// BMA TODO: add convenience functions for repeated matrix ops

/** this function calculates the objective function (negative log
    likelihood) and/or the constraint functions and/or their analytical
    gradients and/or the hessian of the objective function using a 
    precompute and store (store across sequential calls to this function) 
    strategy to reduce the computational cost KRD 2010.05.13
*/
void GradKrigingModel::masterObjectiveAndConstraints(const MtxDbl& theta, int obj_der_mode, int con_der_mode)
{
  // if(obj_der_mode>=1) (1=2^0=> 0th derivative) calculate objective function
  // if(obj_der_mode>=2) (2=2^1=> 1st derivative) calculate objective function and its gradient
  // if(obj_der_mode>=4) (4=2^2=> 2nd derivative) calculate objective function and its gradient and Hessian
  // ERROR if(obj_der_mode>7) ERROR
  //
  // if(con_der_mode>=1) (1=2^0=> 0th derivative) calculate the constraint functions
  // if(con_der_mode>=2) (2=2^1=> 1st derivative) calculate the constraint functions and their gradients
  // ERROR if(con_der_mode>=4) (4=2^2=> 2nd derivative) calculate the constraint functions and their gradients and Hessians
  // ERROR if(con_der_mode>3) ERROR

  //printf("constraintType=|%c| maxConDerMode=%d con_der_mode=%d maxObjDerMode=%d obj_der_mode=%d\n",
  //constraintType,maxConDerMode,con_der_mode,maxObjDerMode,obj_der_mode);

  //might want to replace this with a thrown exception
  assert(//(maxObjDerMode<=7)&&(maxConDerMode<=3)&& ///I've yanked out the analytical derivatives because something in them wasn't working I need to go back and fix this but for now don't let the user use this
	 (maxObjDerMode<=1)&&(maxConDerMode<=1)&&
	 (0<=obj_der_mode)&&(obj_der_mode<=maxObjDerMode)&&
	 (0<=con_der_mode)&&(con_der_mode<=maxConDerMode)&&
	 ((1<=obj_der_mode)||(1<=con_der_mode)));

  //next 5 or 6 must be copied out (can't allow external code to change these or have no guarantee that we can reuse "previous" values)
  //member variable: double obj;
  //member variable: MtxDbl gradObj(numTheta);
  //member variable: MtxDbl hessObj(numTheta,numTheta);
  //member variable: MtxDbl con(numConFunc);
  //member variable: MtxDbl gradCon(numConFunc,numTheta);
  //???member variable: vector<MtxDbl(numTheta,numTheta)> hessCon(numConFunc) ..hessians of the constraints are also a bit more expensive to copy out, might be better to use finite difference hessians of the constraint functions, analytical hessians of the constraint functions are not yet implemented
  //end of variables to copy out

  //these 2 are set by choice of optimizer and are not changed inside this function
  //member variable: int maxObjDerMode;
  //member variable: int maxConDetMode;

  //these are temporary variables that we keep around during emulator 
  //creation so we don't have to constantly allocate and deallocate them, 
  //these get deallocated at the end of create
  //member variable: MtxDbl temp(nTrend)
  //member variable: MtxDbl temp2(numPoints)
  //member variable: MtxDbl temp3(numPoints)
  //member variable: MtxDbl temp4(numPoints)
  //member variable: MtxDbl temp5(numPoints)


  //these are private and their values need to be retained between sequential calls to masterObjectiveAndConstraints, that is no other function (other than the create) can access them, these get deallocated at the end of create
  //member variable: int prevObjDerMode
  //member variable: int prevConDerMode
  //member variable: MtxDbl prevTheta(1,numTheta)
  //member variable: MtxDbl allEigVect(numPoints,numPoints)
  //member variable: MtxDbl allEigVal(numPoints)
  //member variable: MtxDbl Z(numPoints*numPoints,numTheta)
  //member variable: MtxDbl R(numPoints,numPoints)
  //member variable: MtxDbl G(numPoints,nTrend)


  //keep these around after emulator creation so we can evaluate the emulator
  //member variable: MtxDbl rhs(numPoints)
  //member variable: MtxDbl RLU(numPoints,numPoints)
  //member variable: MtxInt ipvt_RLU
  //member variable: MtxDbl Rinv(numPoints,numPoints) //needed to eval integral of adjusted variance, capability not yet added
  //member variable: MtxDbl Rinv_G(numPoints,nTrend)
  //member variable: MtxDbl Gtran_Rinv_G_LU(nTrend,nTrend)
  //member variable: MtxInt ipvt_Gtran_Rinv_G_LU(nTrend)

  //if theta was the same as the last time we called this function than we can reuse some of the things we calculated last time
  
  int i;

  //MtxDbl temp_row_vec_to_print_corrlen(1,numVarsr);
  //for(i=0; i<numTheta; ++i)
  //temp_row_vec_to_print_corrlen(0,i)=1.0/sqrt(2.0*theta(0,i));
  //scaler.unScaleXrDist(temp_row_vec_to_print_corrlen);
  //printf("L=(%g",temp_row_vec_to_print_corrlen(0,0));
  //for(i=1; i<numTheta; ++i)
  //printf(",%g",temp_row_vec_to_print_corrlen(0,i));
  //printf(")\n");


  if(prevTheta.getNElems()!=numTheta) {
    prevTheta.newSize(1,numTheta);
    prevObjDerMode=prevConDerMode=0; 
  }
  else
    for(i=0; i<numTheta; ++i) 
      if(prevTheta(0,i)!=theta(0,i)) {
	prevObjDerMode=prevConDerMode=0;
	break;
      }

  if((obj_der_mode<=prevObjDerMode)&&
     (con_der_mode<=prevConDerMode)) {
    //we've already calculated everything you just asked for so reuse it
    return;
  }

  //record the current theta as the previous theta so we can tell if we
  //should reuse the stuff we calculate this time
  //MtxDbl temp_row_vec_to_print_corrlen(1,numVarsr);
  if((prevObjDerMode==0)&&(prevConDerMode==0)) {
    //printf("Changing correlati0n lengths: "); 
    for(i=0; i<numTheta; ++i) {
      //temp_row_vec_to_print_corrlen(0,i)=1.0/sqrt(2.0*theta(0,i));
      prevTheta(0,i)=theta(0,i); 
    }}
  //scaler.unScaleXrDist(temp_row_vec_to_print_corrlen);
  //printf("L=(%g",temp_row_vec_to_print_corrlen(0,0));
  //for(i=1; i<numTheta; ++i)
  //printf(",%g",temp_row_vec_to_print_corrlen(0,i));
  //printf(")\n");

  int chol_info;

  if((prevObjDerMode==0)&&(prevConDerMode==0)) {
    //R.newSize(numRowsR,numRowsR);
    R.newSize(numEqnAvail,numEqnAvail);
    correlation_matrix(theta); //fills member variable R with exp(Z*theta) 
    //where Z is a member variable calculated from XR, plus derivatives 
    //with respect to XR 
    apply_nugget_build(); //modify R by nug in place
  }

  if((prevObjDerMode==0)&&((1<=obj_der_mode)||(1<=con_der_mode))) {


    //perform LU (replaced with Cholesky) decomposition of R and calculate the determinant of R 
    //http://en.wikipedia.org/wiki/Determinant#Determinant_from_LU_decomposition
    //the difference is that det(L) isn't one for Cholesky, instead det(L)=det(U)

    //this starts with the function values and derivatives at the Anchor points and then adds as many other derivative equations as it can, one at a time while keeping rcondR equal to or greater than the minimum that is allowed (i.e. 1.0/maxCondNum)
    //correlations.copy(theta); //for debug only, bob

    equationSelectingPrecondCholR();
    if(rcondR==0.0) {
      obj=HUGE_VAL;
      con.newSize(numConFunc,1);
      for(int i=0; i<numConFunc; ++i)
	con(i,0)=1.0; //say the constraints are violated
      return;
    }

    //RChol.copy(R);
    //chol_info=0;
    //Chol_fact(RChol,chol_info,rcondR); //preconditioned Cholesky, when Kriging
    //double log_determinant_R, sign_of_det_R;
    //protected_pseudo_inverseR(rcondR, log_determinant_R); //this fills Rinv

    //Rinv.copy(R);
    //pseudo_inverse_sym(Rinv, 1.0/maxCondNum, rcondR, log_determinant_R, sign_of_det_R);
    //rcondR=1.0/(maxCondNum-1.0);

    //is not gradient enhaced R won't need/benefit preconditiong since it has 
    //all 1's on the diagonals
    //printf("Chol\n");
    //assert(chol_info==0);  //here for debug, decide what to do about it later

    double log_determinant_R = 0.0; //need to do this to avoid underflow error
    //for large numbers of points (read as large number of Rows in R), 
    //log(0)=-inf
    for (int i = 0; i < numRowsR; ++i) 
      log_determinant_R += std::log(RChol(i,i)); 
    log_determinant_R *= 2.0; //only multiply by 2 for Cholesky factorization 
    //of R because det(L)=det(U) and det(R)=det(L)*det(U)=det(L)^2
    //so log(det(R))=2*log(det(L))

    //determinant_R=fabs(determinant_R); //KRD added fabs for LU factorization
    //because "The determinant of a positive definite matrix is always positive" http://mathworld.wolfram.com/PositiveDefiniteMatrix.html and det(R)=det(pivot Mtx)*det(L)*det(U); det(L)=1, det(U) is what we calculated above and det(pivot Mtx)=+/- 1, left this comment in, in case someone decides to switch back to LU decomposition

    // precompute, O(M^3) ops, and store quantities so(when give the 
    // derivative of the correlation matrix, R, with respect to a scalar 
    // correlation) we can later evaluate each derivative ok the objective 
    // function (-log(likelihood)) with respect to a scalar component of 
    // correlation using only 
    //   *a Trace(matrix * matrix multiplication) O(M^2) ops
    //   *matrix vector multiplication O(M^2) ops
    //   *vector-vector addition O(M) ops
    //   *dot products O(M) ops
    //   *scalar division, multiplication and addition O(1) ops
    //   
    //  the math
    //  R=(1+nug)^-1*(R+nug*I), typically nug=0
    //  betaHat=(G^T*R^-1*G)^-1*(G^T*R^-1*Y)
    //  eps=Y-G*betaHat
    //  estVarMLE=1/N*eps^T*R^-1*eps
    //  obj=0.5*[N*(log(2*pi)+log(estVarMLE)+1)+log(det(R))]
    //
    //  for gradients of the objective function, precompute and store
    //  G^T*R^-1
    //  (G^T*R^-1*G)^-1
    //  R^-1*G*betaHat-R^-1*Y
    //  R^-1*eps
    //
    //  will loop over k to fill up the gradient of the objective function
    
    //Do the generalized (by R^-1) least squares using min # of ops
    Rinv_G.newSize(numRowsR,nTrend); //precompute and store
    //matrix_mult(Rinv_G,Rinv,G);
    solve_after_Chol_fact(Rinv_G,RChol,G);

  
    //Gtran_Rinv_G_inv.newSize(nTrend,nTrend);
    //matrix_mult(Gtran_Rinv_G_inv,G,Rinv_G,0.0,1.0,'T','N');
    //double log_determinant_Gtran_Rinv_G=0.0, sign_of_det_Gtran_Rinv_G=1.0;
    //pseudo_inverse_sym(Gtran_Rinv_G_inv, 1.0/maxCondNum, rcond_Gtran_Rinv_G, log_determinant_Gtran_Rinv_G, sign_of_det_Gtran_Rinv_G);


    Gtran_Rinv_G_Chol.newSize(nTrend,nTrend);
    matrix_mult(Gtran_Rinv_G_Chol,G,Rinv_G,0.0,1.0,'T','N');
    //double rcond_Gtran_Rinv_G;
    Chol_fact_workspace(Gtran_Rinv_G_Chol,Gtran_Rinv_G_Chol_Scale,Gtran_Rinv_G_Chol_DblWork,Gtran_Rinv_G_Chol_IntWork,chol_info,rcond_Gtran_Rinv_G);
    //printf("rcond_Gtran_Rinv_G = %g\n",rcond_Gtran_Rinv_G);

    //if(~(chol_info==0)) assert(chol_info==0);  //for debug, do something else for production




    double log_determinant_Gtran_Rinv_G=0.0; //need this for the Rassmussen & Williams formulation of likelihood
    for (int itrend = 0; itrend < nTrend; ++itrend)
      log_determinant_Gtran_Rinv_G += std::log(Gtran_Rinv_G_Chol(itrend,itrend)); 
    log_determinant_Gtran_Rinv_G *= 2.0; //only for Cholesky factorization of R

    temp.newSize(nTrend,1);
    matrix_mult(temp, Rinv_G, Y, 0.0, 1.0, 'T', 'N');
    betaHat.newSize(nTrend,1);
    solve_after_Chol_fact(betaHat,Gtran_Rinv_G_Chol,temp); //O(nTrend^2) ops
    //matrix_mult(betaHat,Gtran_Rinv_G_inv,temp); //O(nTrend^2) ops

    temp2.copy(Y); //this will be eps=epsilon=Y-G(XR)*betaHat, but use 
    //variable temp2 because we would only need the variable "eps" for these 
    //5 lines of code (not counting comments) and we want to save space, 
    //afterwards we will only need R^-1*eps which is stored in "rhs"
    matrix_mult(temp2, G, betaHat, 1.0, -1.0, 'N', 'N'); //eps=Y-G(XR)*betaHat
    rhs.newSize(numRowsR,1);
    //matrix_mult(rhs,Rinv,temp2);
    solve_after_Chol_fact(rhs,RChol,temp2);

    /*
    MtxDbl GbetaHat;
    matrix_mult(GbetaHat,G,betaHat);
    MtxDbl y(GbetaHat);
    matrix_mult(y,R,rhs,1.0,1.0);

    FILE *fp=fopen("test_scaled_Y2.txt","w");
    fprintf(fp,"XR=");
    for(int i=0; i<XR.getNRows(); ++i) {
      fprintf(fp,"\n%-12.6g",XR(i,0));
      for(int j=1; j<XR.getNCols(); ++j)
	fprintf(fp," %-12.6g",XR(i,j));
    }
    
    fprintf(fp,"\n\nG=");
    for(int i=0; i<numPoints; ++i) {
      fprintf(fp,"\n%-12.6g",G(i,0));
      for(int j=1; j<G.getNCols(); ++j)
	fprintf(fp," %-12.6g",G(i,j));
    }
    
    fprintf(fp,"\n\nbetaHat=");
    for(int i=0; i<betaHat.getNRows(); ++i)
      fprintf(fp,"\n%-12.6g",betaHat(i,0));

    fprintf(fp,"\n\nGbetaHat=");
    for(int i=0; i<numPoints; ++i)
      fprintf(fp,"\n%-12.6g",GbetaHat(i,0));

    fprintf(fp,"\n\nR=");
    for(int i=0; i<numPoints; ++i) {
      fprintf(fp,"\n%-12.6g",R(i,0));
      for(int j=1; j<R.getNCols(); ++j)
	fprintf(fp," %-12.6g",R(i,j));
    }
  
    fprintf(fp,"\n\nrhs=");
    for(int i=0; i<rhs.getNRows(); ++i)
      fprintf(fp,"\n%-12.6g",rhs(i,0));

    fprintf(fp,"\n\n[scaled_Y y]=");
    for(int i=0; i<numPoints; ++i)
      fprintf(fp,"\n%-12.6g %-12.6g",Y(i,0),y(i,0));
    fflush(fp);
    fclose(fp);

    fp=fopen("test_emulator_reproduce_y.txt","w");
    int iii;
    fprintf(fp,"R=");
    for(iii=0; iii<R.getNRows(); ++iii) {
      fprintf(fp,"\n%-12.6g",R(iii,0));
      for(int jjj=1; jjj<R.getNCols(); ++jjj)
	fprintf(fp," %-12.6g",R(iii,jjj));
    }

    //#ifdef __NO_NO_NO_NO__
    MtxDbl RRR(R);

    fprintf(fp,"\n\nRRR=");
    for(iii=0; iii<RRR.getNRows(); ++iii) {
      fprintf(fp,"\n%-12.6g",RRR(iii,0));
      for(int jjj=1; jjj<RRR.getNCols(); ++jjj)
	fprintf(fp," %-12.6g",RRR(iii,jjj));
    }
	      
    MtxDbl GGG(G);
    MtxDbl YYY(Y);
    int infoRRR, info_GGGtran_RRRinv_GGG;
    double rcondRRR, rcond_GGGtran_RRRinv_GGG;

    fprintf(fp,"\n\nRChol (rcondR=%g)=",rcondR);
    for(iii=0; iii<RChol.getNRows(); ++iii) {
      fprintf(fp,"\n%-12.6g",RChol(iii,0));
      for(int jjj=1; jjj<RChol.getNCols(); ++jjj)
	fprintf(fp," %-12.6g",RChol(iii,jjj));
    }


    Chol_fact(RRR,infoRRR,rcondRRR);

    fprintf(fp,"\n\nRRRChol (rcondRRR=%g)=",rcondRRR);
    for(iii=0; iii<RRR.getNRows(); ++iii) {
      fprintf(fp,"\n%-12.6g",RRR(iii,0));
      for(int jjj=1; jjj<RRR.getNCols(); ++jjj)
	fprintf(fp," %-12.6g",RRR(iii,jjj));
    }

    MtxDbl RRRinv_GGG;
    solve_after_Chol_fact(RRRinv_GGG,RRR,GGG);

    fprintf(fp,"\n\nRinv_G=");
    for(iii=0; iii<Rinv_G.getNRows(); ++iii) {
      fprintf(fp,"\n%-12.6g",Rinv_G(iii,0));
      for(int jjj=1; jjj<Rinv_G.getNCols(); ++jjj)
	fprintf(fp," %-12.6g",Rinv_G(iii,jjj));
    }

    fprintf(fp,"\n\nRRRinv_GGG=");
    for(iii=0; iii<RRRinv_GGG.getNRows(); ++iii) {
      fprintf(fp,"\n%-12.6g",RRRinv_GGG(iii,0));
      for(int jjj=1; jjj<RRRinv_GGG.getNCols(); ++jjj)
	fprintf(fp," %-12.6g",RRRinv_GGG(iii,jjj));
    }


    MtxDbl GGGtran_RRRinv_GGG_Chol;
    matrix_mult(GGGtran_RRRinv_GGG_Chol,GGG,RRRinv_GGG,0.0,1.0,'T','N');
    Chol_fact(GGGtran_RRRinv_GGG_Chol,info_GGGtran_RRRinv_GGG,rcond_GGGtran_RRRinv_GGG);
    MtxDbl TTTemp;
    matrix_mult(TTTemp,RRRinv_GGG,YYY,0.0,1.0,'T','N');
    MtxDbl BBBeta;
    solve_after_Chol_fact(BBBeta,GGGtran_RRRinv_GGG_Chol,TTTemp);
    MtxDbl TTTemp2(Y);
    matrix_mult(TTTemp2,GGG,BBBeta,1.0,-1.0,'N','N');
    MtxDbl rrrhs;
    solve_after_Chol_fact(rrrhs,RRR,TTTemp2);
    correlations.copy(theta);
    MtxDbl rrr;
    correlation_matrix(rrr,XR);
    MtxDbl IIIchop;
    solve_after_Chol_fact(IIIchop,RRR,rrr,'T');
    MtxDbl TTTemp3;
    matrix_mult(TTTemp3,GGG,BBBeta,0.0,1.0,'N','N');
    MtxDbl TTTemp4(TTTemp3);
    matrix_mult(TTTemp3,rrr,rrrhs,1.0,1.0,'N','N');
    matrix_mult(TTTemp4,IIIchop,TTTemp2,1.0,1.0,'T','N');
    MtxDbl ggg;
    eval_trend_fn(ggg,XR);
    MtxDbl TTTemp5;
    matrix_mult(TTTemp5,ggg,BBBeta);
    matrix_mult(TTTemp5,rrr,rrrhs,1.0,1.0,'N','N');
    
    

    fprintf(fp,"numPoints=%d\nnumEqnKeep=%d\nnumRowsR=%d\n\n[Y           YYY          E(YYY3)      E(YYY4)      E(YYY5)]=\n",numPoints,numEqnKeep,numRowsR);
    for(iii=0; iii<TTTemp3.getNRows(); ++iii)
      fprintf(fp,"%-12.6g %-12.6g %-12.6g %-12.6g %-12.6g\n",Y(iii,0),YYY(iii,0),TTTemp3(iii),TTTemp4(iii),TTTemp5(iii));
    for(; iii<YYY.getNRows(); ++iii)
      fprintf(fp,"%-12.6g %-12.6g\n",Y(iii,0),YYY(iii,0));

    fprintf(fp,"\n[betaHat     BBBeta]=\n");
    for(iii=0; iii<nTrend; ++iii)
      fprintf(fp,"%-12.6g %-12.6g\n",betaHat(iii,0),BBBeta(iii,0));


    fprintf(fp,"\n[rhs         rrrhs]=\n");
    for(iii=0; iii<rhs.getNRows(); ++iii)
      fprintf(fp,"%-12.6g %-12.6g\n",rhs(iii,0),rrrhs(iii,0));

    fflush(fp);
    fclose(fp);
    //    assert(false);
    */

    //it's actually the log likelihood, which we want to maximize
    //likelihood = -0.5*(numRowsR*(std::log(4.0*std::acos(0.0))+std::log(estVarianceMLE)+1)
    //		       +std::log(determinant_R)); //from Koehler and Owen 

#ifdef __NKM_UNBIASED_LIKE__
    //derived following: C. E. Rasmussen & C. K. I. Williams, Gaussian Processes for Machine Learning, the MIT Press, 2006, ISBN 026218253X. c 2006 Massachusetts Institute of Technology. www.GaussianProcess.org/gpml...  we assume a "vague prior" (i.e. that we don't know anything) for betaHat, then like "Koehler and Owen" we replace the covariance matrix K with (unadjusted variance)*R (where R is the correlation matrix) and find unadjusted variance and betaHat through maximum likelihood.

    //the unbiased estimate of unadjusted variance
    estVarianceMLE = dot_product(temp2,rhs)/(numRowsR-nTrend); 

    //the "per point" unbiased log(likelihood)
    likelihood = -0.5*(std::log(estVarianceMLE)+(log_determinant_R+log_determinant_Gtran_Rinv_G)/(numRowsR-nTrend)); 

#else
    //derived the "Koehler and Owen" way (assumes we know the trend function, and is therefore biased, but usally seems to work better for surrogate based optimization)

    //the estimate of unadjusted variance
    estVarianceMLE = dot_product(temp2,rhs)/numRowsR; //the "Koehler and Owen" way

    //the "per point" log(likelihood), with constant terms dropped off
    likelihood = -0.5*(std::log(estVarianceMLE)+log_determinant_R/numRowsR);
    
    //"per point" log(likelihood) with constant terms dropped off times the number of points... poor performance (and reproducing true surfaces)
    //likelihood = -0.5*(numRowsR*std::log(estVarianceMLE)+log_determinant_R); 

    //likelihood = -0.5*(numRowsR*(std::log(4.0*std::acos(0.0))+std::log(estVarianceMLE)+1.0)
    //  		       +log_determinant_R); //from Koehler and Owen 

#endif

    //if(likelihood>=DBL_MAX)
    //printf("[estVarianceMLE=%g determinant_R=%g]",estVarianceMLE,determinant_R);

    //the objective function being MINIMIZED is the negative of the log 
    //likelihood (on a per point basis so numbers will be comparable 
    //regardless of how many points there are)
    obj=-likelihood;  
    //obj=likelihood;  
    //printf("[obj=%g]",obj);

    prevObjDerMode=1; //increase prevObjDerMode to the current value
    if((obj_der_mode==1)&&(con_der_mode<=prevConDerMode)) {
      //we have everything we need so exit early
      return;
    }
  }


  if((prevConDerMode==0)&&(1<=con_der_mode)) {
    //calculate the constraint that ensures that the correlation matrix is 
    //well conditioned.  conceputally: condition_number < number_of_points; 
    //implemented as: (largest_eigenvalue+nug)/ number_of_points - (smallest_eigenvalue(s)+nug) < 0 
    //(so we don't have to worry about dividing by zero when the smallest eigenvalue is zero 
    //and the scaling should be rougly consistant regardless of the number of points), the 
    //value of nug is typically zero
    con.newSize(numConFunc,1);
    
    if(constraintType.compare("rcond")==0) { //use rcond (and maybe its numerical derivatives) to bound the 
      //condition number
      
      assert((1<=prevObjDerMode)&&(numConFunc==1)); //make sure we have calculated rcondR already
      //con(0,0)=1.0-rcondR*3.0*maxCondNum;  //have seen rcond as low as about 1/3 of the true value
      con(0,0)=1.0-rcondR*maxCondNum;  //have seen rcond as low as about 1/3 of the true value
    }
    /*
    else if(constraintType.compare("eig")==0) { //use eigenvalues (and maybe their analytical derivatives) to
      //bound the condition number
      allEigVect.newSize(numPoints,numPoints);
      allEigVal.newSize(numPoints);

      eig_sym(allEigVect, allEigVal, R); //R is correct
      for(int icon=0; icon<numConFunc; icon++)
	con(icon,0)=(allEigVal(numPoints-1,0)+nug)/maxCondNum-(allEigVal(icon,0)+nug);
        //con(icon,0)=1.0-maxCondNum*(allEigVal(icon,0)+nug)/(allEigVal(numPoints-1,0)+nug);
        //con(icon,0)=1.0/maxCondNum-(allEigVal(icon,0)+nug);
    }
    else
      assert((constraintType.compare("eig")==0)||(constraintType.compare("rcond")==0));
    */
    else
      assert(constraintType.compare("rcond")==0);
    

    prevConDerMode=1; //increase prevConDerMode to current value
    if((con_der_mode==1)&&(obj_der_mode<=prevObjDerMode)) {
      //we have everything we need so exit early
      return;
    }
  }
  
  
  //******************************************************************
  //******************************************************************
  // calculating the gradients of the constraint functions (when we 
  // don't need the gradient of the objective function) starts here
  //******************************************************************
  //******************************************************************
  if(((obj_der_mode<2)||(2<=prevObjDerMode))&&
     ((prevConDerMode<2)&&(2<=con_der_mode))) {
    assert(false);
  }

  //******************************************************************
  //******************************************************************
  // calculating the gradient of the objective function (and optionally
  // the constraint functions) starts here
  //******************************************************************
  //******************************************************************
  if((prevObjDerMode<2)&&(2<=obj_der_mode)) {
    assert(false);
  }
  

  //******************************************************************
  //******************************************************************
  // calculating the hessian of the objective function starts here
  //******************************************************************
  //******************************************************************
  if((prevObjDerMode<4)&&(4<=obj_der_mode)) {
    assert(false);
  }//if((prevObjDerMode<4)&&(4<=obj_der_mode))

  return;
}


void GradKrigingModel::getRandGuess(MtxDbl& guess) const
{
  int mymod = 1048576; //2^20 instead of 10^6 to be kind to the computer
  guess.newSize(1,numVarsr);
  //double corr_length,tempdouble;
  for(int j=0;j<numVarsr;j++) {
    //tempdouble=(std::rand() % mymod);
    //tempdouble/=mymod;
    //printf("tempdouble=%g ",tempdouble);
    //corr_length=std::pow(2.0,tempdouble*(log2(max_corr_length)-log2(min_corr_length)) + log2(min_corr_length));    
      
    //corr_length=std::exp((std::rand() % mymod)*(maxNatLogCorrLen-minNatLogCorrLen)/mymod+
    //minNatLogCorrLen);
    //guess(j) = 1.0/(2.0*corr_length*corr_length); 
    guess(0,j) = (std::rand() % mymod)*(maxNatLogCorrLen-minNatLogCorrLen)/mymod+
      minNatLogCorrLen; //this returns a random nat_log_corr_len which is the space we need to search in

    //guess(0,j)=100.0;
  }
  //printf("\n");
}


/** this functions makes guessed values of the correlation paramters
    feasible, i.e. decreases the condition number of the correlation
    matrix to be less than the number of points by taking steps in the
    direction opposite the gradient of
    objective_function=largest_eigenvalue_of_R/number_of_points -
    smallest_eigenvalue_of_R with a step size such that the
    objective_function would be reduced to zero in a single step if it
    were a linear function, R really is R in this function, it determines
    the needed/desired nugget so R can be modified to R OUTSIDE of this 
    function */
MtxDbl& GradKrigingModel::makeGuessFeasible(MtxDbl& nat_log_corr_len, OptimizationProblem *opt) {
  int k;
  //int chol_info;
  MtxDbl theta(1,numTheta);
  for(k=0; k<numTheta; ++k)
    theta(0,k)=0.5*std::exp(-2.0*nat_log_corr_len(0,k));

  //R.newSize(numRowsR,numRowsR);
  R.newSize(numEqnAvail,numEqnAvail);
  
  //temp.newSize(numRowsR);
  temp.newSize(numEqnAvail,1);
  correlation_matrix(theta); //assigns to member variable R

  if((ifChooseNug==true)||(nug<0.0))
    nug=0.0;

  apply_nugget_build();

  //RChol.copy(R);
  double best_rcond; //, log_determinant_R;
  equationSelectingPrecondCholR();
  best_rcond=rcondR;
  //Chol_fact(RChol,chol_info,best_rcond);
  //protected_pseudo_inverseR(best_rcond, log_determinant_R);
  //if(~(chol_info==0)) assert(chol_info==0); //for debug, do something different for production

  double rcond;
  MtxDbl guess(1,numTheta);
  MtxDbl guess_theta(1,numTheta);
  int iguess=0;
  while((1.0-best_rcond*maxCondNum>0.0)&&(iguess<50)) {
    iguess++;
    getRandGuess(guess);

    //convert guess from nat_log_corr_len to theta
    for(k=0; k<numTheta; ++k)
      guess_theta(0,k)=0.5*std::exp(-2.0*guess(0,k));

    correlation_matrix(guess_theta); //assigns to member variable R
    apply_nugget_build();
    //RChol.copy(R);
    //Chol_fact(RChol,chol_info,rcond);
    //protected_pseudo_inverseR(rcond, log_determinant_R);
    equationSelectingPrecondCholR();
    rcond=rcondR;

    //if(~(chol_info==0)) assert(chol_info==0); //for debug, do something else for production

    if(rcond>best_rcond) {
      
      best_rcond=rcond;
      theta.copy(guess_theta);
      nat_log_corr_len.copy(guess);
    }
  }
  
  if(constraintType.compare("rcond")!=0) 
    assert(false);
    
  return nat_log_corr_len; //we can't compute analytical derivatives of rcond only the eigenvalues
}

/** evaluate the trend function g(x), using specified {Poly, Rot}
Here, g() and x can represent arbitrary (elsewhere represented by
g(x)) or data points (elsewhere represented by G(X)).  The trend
function (unadjusted mean) is dot(g(x),betHat), this returns (matrix)
g(x) for collection of points x. KRD originally implemented this with
a linear trend function */
MtxDbl& GradKrigingModel::
eval_trend_fn(MtxDbl& g, const MtxInt& poly, const MtxInt& der,  
	      const MtxDbl& rot_or_eul_ang, const MtxDbl& x) const
{
  //printf(" you called me ");
  MtxDbl xx;
  rotate_xr(xx, rot_or_eul_ang, x);
  //int numVars=x.getNCols(); int npts=x.getNRows();
  //for(int i=0; i<npts; i++) {
  //printf("x(%d,:)=[ ",i);
  //for(int j=0; j<numVars; j++)
  //printf("%8f ",x(i,j));
  //printf("];  xx(%d,:)=[ ",i);
  //for(int j=0; j<numVars; j++)
  //printf("%8f ",xx(i,j));
  //printf("]\n");
  //}
  return (evaluate_poly_der_basis(g, poly, der, xx));
  ///return (LinearRegressionModel::evalBasis(g,poly,rot_or_eul_ang,x));
}


/// evaluate the trend function g(xr), using class members {Poly, Rot} and only the zeroth order derivative
MtxDbl& GradKrigingModel::eval_trend_fn(MtxDbl& g, const MtxDbl& xr) const
{
  //MtxDbl& evalTrendFunc(MtxDbl& g, MtxDbl& xr) {
  //printf("KMeTF <");
  
  MtxInt der(1,numVarsr); der.zero();
  eval_trend_fn(g, Poly, der, Rot, xr);
  //printf("<KMeTF\n");
  return (g);
    //return (LinearRegressionModel::evalBasis(g,Poly,Rot,xr));
}

// BMA TODO: These need to be moved to optimizer and then any defauls
// overridden here

void GradKrigingModel::set_conmin_parameters(OptimizationProblem& opt) const
{
  //set conmin specific parameters for this problem
  //in dot 4.2 the analytical derivative order of objective and constraints must be the same

  assert((maxObjDerMode==1)&&(maxConDerMode==1)); //no analytical gradients

  if((maxObjDerMode==1)&&(maxConDerMode==1))
    opt.conminData.nfdg = 0; //use numerical  gradients of objective and constraints
  else if((maxObjDerMode==3)&&(maxConDerMode==3))
    opt.conminData.nfdg = 1; //use analytical gradients of objective and constraints 
  else if((maxObjDerMode==3)&&(maxConDerMode==1))
    opt.conminData.nfdg = 2; //uses analytical derivatives for the objective and numerical derivatives for constraints
  else
    assert(false);

  opt.conminData.iprint = 0; //ammount of to screen output from Conmin
  opt.conminData.itmax  = maxTrials; //maximum # of Conmin iterations
  opt.conminData.fdch   = 1.0e-2; //Relative finite difference step size.
  opt.conminData.fdchm  = 1.0e-2; //Absolute finite difference step size.
  opt.conminData.ct     = -0.1; // Constraint thickness parameter, The absolute value of CT decreases in magnitude during optimization. 
  opt.conminData.ctmin  = 0.004; //Minimum absolute value of CT used during optimization.
  opt.conminData.ctl    = -0.01; //Constraint thickness parameter for linear and side constraints.
  opt.conminData.ctlmin = 0.001; //Minimum value of CTL used during optimization.
  opt.conminData.delfun = 0.001; //Relative convergence criterion threshold, Threshold for the minimum relative change in the objective function
  opt.conminData.dabfun = 0.001; //Absolute convergence criterion threshold. Threshold for the minimum relative change in the objective function
  opt.conminData.nside  = 1; //side constraints parameter
  opt.conminData.itrm   = 3; //diminishing return criterion iteration number
  opt.conminData.icndir = numTheta+1; //conjugate direction restart parameter
}

void GradKrigingModel::set_direct_parameters(OptimizationProblem& opt) const
{
  opt.directData.minBoxSize = -1.0;
  opt.directData.volBoxSize = -1.0;
  //opt.directData.minBoxSize = 1.0e-15;
  //opt.directData.volBoxSize = 1.0e-15;
  //opt.directData.minBoxSize = 1.0e-3;
  //opt.directData.volBoxSize = 1.0e-5;
  opt.directData.solutionTarget = -DBL_MAX;
  opt.directData.convergenceTol = 1.0e-4;
  opt.directData.maxFunctionEvals = maxTrials;
  opt.directData.maxIterations = 1000; 
  opt.directData.verboseOutput = false;
  opt.directData.constraintsPresent = true;
}

} // end namespace nkm
