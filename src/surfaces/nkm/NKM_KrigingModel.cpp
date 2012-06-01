#include "NKM_SurfPack.hpp"
#include "NKM_KrigingModel.hpp"
//#include "Accel.hpp"
#include "NKM_LinearRegressionModel.hpp"
#include <math.h>
#include <iostream>
#include <cfloat>

namespace nkm {

using std::cout;
using std::cerr;
using std::endl;
using std::ostringstream;

//#define __NKM_UNBIASED_LIKE__

//#define __KRIGING_DER_TEST__


void KrigingModel::preAllocateMaxMemory() {
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

void KrigingModel::nuggetSelectingCholR(){
  numEqnKeep=numRowsR=numPoints;
  iEqnKeep.newSize(numEqnAvail,1);
  for(int i=0; i<numEqnKeep; ++i)
    iEqnKeep(i,0)=i;
  nug=0.0;
  polyOrder=polyOrderRequested;
  while((numRowsR<=numTrend(polyOrder,0))&&(polyOrder>0))
    --polyOrder;
  nTrend=numTrend(polyOrder,0);
  for(int i=0; i<numEqnKeep; ++i) 
    Y(i,0)=Yall(iEqnKeep(i,0),0);

  for(int itrend=0; itrend<nTrend; ++itrend) 
    for(int i=0; i<numEqnKeep; ++i) 
      G(i,itrend)=Gall(iEqnKeep(i,0),itrend);

  double min_allowed_rcond=1.0/maxCondNum;
  RChol.copy(R);
  int ld_RChol=RChol.getNRowsAct();
  int chol_info;
  scaleRChol.newSize(numEqnAvail,1); //no scaling is necessary since have
  //all ones on the diagonal of R but the generic preconditioned cholesky 
  //factoriztion function doesn't know that in advance
  rcondDblWork.newSize(3*ld_RChol,1);
  rcondIntWork.newSize(ld_RChol,1);
  Chol_fact_workspace(RChol,scaleRChol,rcondDblWork,rcondIntWork,
		      chol_info,rcondR);
  if(rcondR<=min_allowed_rcond) {
    double dbl_num_pts=static_cast<double>(numPoints);
    double sqrt_num_pts=std::sqrt(dbl_num_pts);
    min_allowed_rcond*=sqrt_num_pts;
    rcondR/=sqrt_num_pts; //one norm is within a factor of N^0.5 of 2 norm
    double min_eig_worst=(rcondR*dbl_num_pts)/(1.0+(dbl_num_pts-1.0)*rcondR);
    double max_eig_worst=dbl_num_pts-(dbl_num_pts-1.0)*min_eig_worst;
    nug=(min_allowed_rcond*max_eig_worst-min_eig_worst)/
      (1.0-min_allowed_rcond);
    //this nugget will make the worst case scenario meet (with an ==) the maxCondNum constraint, I (KRD) don't expect this to ever == fail because I don't expect rcond to be *N^-0.5 without nugget and be *N^0.5 with nugget while the maximum eigen value of R (without nugget) is N-(N-1)*min_eigval (all eigenvalues except the largest are the smallest possible for the given rcond) note that rcond is the LAPACK ESTIMATE of the 1 norm condition number so there are no 100% guarantees.
    apply_nugget_build();
    RChol.copy(R);

    Chol_fact_workspace(RChol,scaleRChol,rcondDblWork,rcondIntWork,
			chol_info,rcondR);
  }
  return;
}


void KrigingModel::equationSelectingCholR(){ 
  polyOrder=polyOrderRequested;
  
  //printf("Entered equationSelectingCholR()\n");
  double min_allowed_rcond=1.0/maxCondNum;
  //printf("min_allowed_rcond=%g\n",min_allowed_rcond);
  //exit(0);
  //double min_allowed_pivot_est_rcond=256.0/maxCondNum;


  RChol.copy(R);  //This is of size numEqnAvail by numEqnAvail
  int ld_RChol=RChol.getNRowsAct();
  int chol_info;
  scaleRChol.newSize(numEqnAvail,1); //no scaling is necessary since have
  //all ones on the diagonal of R but the generic preconditioned cholesky 
  //factoriztion function doesn't know that in advance
  rcondDblWork.newSize(3*ld_RChol,1);
  rcondIntWork.newSize(ld_RChol,1);
  Chol_fact_workspace(RChol,scaleRChol,rcondDblWork,rcondIntWork,
		      chol_info,rcondR);
  //printf("rcondR=%g",rcondR);
  iEqnKeep.newSize(numEqnAvail,1);
  if(min_allowed_rcond<rcondR) {
    //printf("numEqnAvail=%d\n",numEqnAvail);
    numRowsR=numEqnKeep=numEqnAvail;
    G.copy(Gall);
    Y.copy(Yall);
    for(int i=0; i<numPoints; ++i)
      iEqnKeep(i,0)=i;
    //printf("\n");
    while((numRowsR<=numTrend(polyOrder,0))&&(polyOrder>0))
      --polyOrder;
    nTrend=numTrend(polyOrder,0);
    return;
  }

  RChol.copy(R);
  if(ifHaveAnchorPoint&&(iAnchorPoint!=0)) {
    //printf("iAnchorPoint=%d\n",iAnchorPoint);
    double dswap;
    for(int i=0; i<numPoints; ++i) {
      dswap=RChol(i,0);
      RChol(i,0)=RChol(i,iAnchorPoint);
      RChol(i,iAnchorPoint)=dswap;
    }
    for(int j=0; j<numPoints; ++j) {
      dswap=RChol(0,j);
      RChol(0,j)=RChol(iAnchorPoint,j);
      RChol(iAnchorPoint,j)=dswap;
    }    
  }

  ld_RChol=RChol.getNRowsAct();
  rcondDblWork.newSize(3*ld_RChol,1);
  rcondIntWork.newSize(ld_RChol,1);

  //if the user specifies an anchor point it must be the first point
  int info=0;
  char uplo='B'; //'B' means we have both halves of R in RChol so the 
  //fortran doesn't have to copy one half to the other
  
  numEqnKeep=numEqnAvail; //=numPoints
  //printf("b4PC numEqnKeep=%d",numEqnKeep);
  PIVOTCHOL_F77(&uplo, &numEqnAvail, RChol.ptr(0,0), &ld_RChol,
    		iEqnKeep.ptr(0,0), &numEqnKeep, &min_allowed_rcond, 
		//&min_allowed_pivot_est_rcond, 
		&info); 
  //printf("Pivoting Cholesky: says numEqnAvail=%d numEqnKeep=%d info=%d",numEqnAvail,numEqnKeep,info);
  //exit(0);
  //if(numEqnKeep<numNeededEqn) {
  //numEqnKeep=numNeededEqn;

  //for(int i=0; i<numEqnAvail; ++i) 
  //iEqnKeep(i,0)-=1;
  
  if(ifHaveAnchorPoint&&(iAnchorPoint!=0)) {
    //printf("have nonzero AnchorPoint\n");
    //convert from fortran indices (starts from 1) to C indices (starts from 0)
    //BUT we swapped point zero and point iAnchorPoint before doing the pivoted Cholesky
    //we have to make sure that change is reflected in iEqnKeep

    iEqnKeep(0,0)=iAnchorPoint; //since there are all ones on the diagonal the first 
    //element would not be pivoted, so we know that this is the AnchorPoint we swapped 
    //to the first/zeroth spot before doing the pivoted Cholesky

    for(int i=1; i<numPoints; ++i) {      
      iEqnKeep(i,0)-=1; 
      if(iEqnKeep(i,0)==iAnchorPoint) //the point that was in AnchorPoint position is
	//really the first/zeroth point
	iEqnKeep(i,0)=0;      
    }
  }
  else {
    //convert from fortran indices (starts from 1) to C indices (starts from 0)    
    for(int i=0; i<numPoints; ++i) 
      iEqnKeep(i,0)-=1;
  }
  
  oneNormR.newSize(numEqnAvail,1);
  sumAbsColR.newSize(numEqnAvail,1); //used in computing the one norm 

  //precompute and store the one norm after adding every equation 
  //in the pivoted Cholesky order, this is prepwork for the 
  //bisection search
  int jeqn=iEqnKeep(0,0);
  for(int i=0; i<numEqnAvail; ++i) {
    sumAbsColR(i,0)=std::fabs(R(iEqnKeep(i,0),jeqn));
  }
  oneNormR(0,0)=sumAbsColR(0,0);
    
  double tempdouble;
  for(int j=1; j<numEqnKeep; ++j) {
    jeqn=iEqnKeep(j,0);
    for(int i=0; i<numEqnAvail; ++i) 	
      sumAbsColR(i,0)+=std::fabs(R(iEqnKeep(i,0),jeqn));
    tempdouble=sumAbsColR(0,0);
    for(int i=1; i<=j; ++i) 
      if(tempdouble<sumAbsColR(i,0))
	tempdouble=sumAbsColR(i,0);
    oneNormR(j,0)=tempdouble;
    //printf("oneNormR(%d,0)=%g\n",j,oneNormR(j,0));
  }
  //exit(0);
  //at this point numEqnKeep is not YET the final number of equations 
  //we will keep, only the maximum number of candidates to keep
  //
  //we will use the LAPACK rcond function with a bisection search to
  //find the last equation that isn't too poorly conditioned, but 
  //first we need the rcond at the point the Cholesky decomposition 
  //terminated.

  uplo='L';
  rcondDblWork.newSize(3*ld_RChol,1);
  rcondIntWork.newSize(ld_RChol,1);
  int icurr_lapack_rcondR=numEqnKeep-1;
  DPOCON_F77(&uplo,&numEqnKeep,RChol.ptr(0,0),&ld_RChol,
	     oneNormR.ptr(icurr_lapack_rcondR,0),
	     &rcondR,rcondDblWork.ptr(0,0),rcondIntWork.ptr(0,0),
	     &info);  
  lapackRcondR.newSize(numEqnAvail,1);
  lapackRcondR.zero(); 
  lapackRcondR(0,0)=1.0; 
  lapackRcondR(icurr_lapack_rcondR,0)=rcondR;
  

  int inext_lapack_rcondR=icurr_lapack_rcondR;
  int iprev_lapack_rcondR=1;
  if(rcondR<=min_allowed_rcond) {
    icurr_lapack_rcondR=numTrend(polyOrder,0);
    int num_needed_eqn=icurr_lapack_rcondR+1;
    DPOCON_F77(&uplo,&num_needed_eqn,RChol.ptr(0,0),&ld_RChol,
	       oneNormR.ptr(icurr_lapack_rcondR,0),
	       &rcondR,rcondDblWork.ptr(0,0),rcondIntWork.ptr(0,0),
	       &info);  
    lapackRcondR(icurr_lapack_rcondR,0)=rcondR;

    if((rcondR==min_allowed_rcond)||
       ((min_allowed_rcond<rcondR)&&(inext_lapack_rcondR==iprev_lapack_rcondR+1))
       ) {
      numEqnKeep=num_needed_eqn; //number of trend function +1 set above
    }
    else{      
      if(rcondR<min_allowed_rcond) {
	inext_lapack_rcondR=icurr_lapack_rcondR;
	iprev_lapack_rcondR=1;
      }
      else{
	iprev_lapack_rcondR=icurr_lapack_rcondR;
      }

      //need at most ceil(log2()) more calls to rcond... this many
      int max_rcond_iter=
	std::ceil(std::log(static_cast<double>
			   (inext_lapack_rcondR-iprev_lapack_rcondR))
		  /std::log(2.0));
        
      //do the bisection search if necessary, at most ceil(log2()) more 
      //calls to rcond
      int rcond_iter=0;
      while((lapackRcondR(inext_lapack_rcondR,0)<=min_allowed_rcond)&&
	    (inext_lapack_rcondR>iprev_lapack_rcondR)) {
	++rcond_iter;
	icurr_lapack_rcondR=(iprev_lapack_rcondR+inext_lapack_rcondR)/2;
	numEqnKeep=icurr_lapack_rcondR+1;
	
	DPOCON_F77(&uplo,&numEqnKeep,RChol.ptr(0,0),&ld_RChol,
		   oneNormR.ptr(icurr_lapack_rcondR,0),
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
	  numEqnKeep=icurr_lapack_rcondR+1;
	  break;
	}
	if((inext_lapack_rcondR-iprev_lapack_rcondR==1)||
	   (max_rcond_iter<rcond_iter)) {
	  numEqnKeep=iprev_lapack_rcondR+1;
	  rcondR=lapackRcondR(iprev_lapack_rcondR,0);
	  break;
	}
      }
    }
  }
  //printf(" pivoted_rcondR=%g numRowsR=%d\n",rcondR,numEqnKeep);
  
  numRowsR=numEqnKeep;
  polyOrder=polyOrderRequested;
  while((numRowsR<=numTrend(polyOrder,0))&&(polyOrder>0))
    --polyOrder;
  nTrend=numTrend(polyOrder,0);
  //printf("nTrend=%d numRowsR=%d\n",nTrend,numRowsR);
  Poly.resize(nTrend,numVarsr); //I am relying on the matrix class's actual 
  //size not changing and that it's contents aren't being written over, so 
  //that when I enlarge the matrix up to polyOrderRequested I recover all 
  //the polynomial powers up to polyOrderRequested

  RChol.resize(numRowsR,numRowsR); //resize() instead of newSize() 
  //because we want to keep the current contents in the same 2D 
  //order
  
  Y.newSize(numEqnKeep,1); //newSize() because we don't care about 
  //the current contents of Y
  G.newSize(numEqnKeep,nTrend);
  iEqnKeep.newSize(numEqnKeep,1);

  for(int i=0; i<numEqnKeep; ++i) 
    Y(i,0)=Yall(iEqnKeep(i,0),0);

  for(int itrend=0; itrend<nTrend; ++itrend) 
    for(int i=0; i<numEqnKeep; ++i) 
      G(i,itrend)=Gall(iEqnKeep(i,0),itrend);
    
}


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
KrigingModel::KrigingModel(const SurfData& sd, const ParamMap& params)
  : SurfPackModel(sd,sd.getJOut()), numVarsr(sd.getNVarsr()), 
    numTheta(numVarsr), numPoints(sd.getNPts()), XR(sdBuild.xr), 
    Yall(sdBuild.y)
{
  //Yall.copy(Y);
  numEqnAvail=numPoints;

  //printf("calling the right KrigingModel constructor\n"); fflush(stdout);

  //if the SurfDataScaler class does what it's supposed to (the only private content in sdBuild that a Model can access are the scaling bits, and only the model can see inside the scaler) the next line will cause an error when you try to compile with it uncommented
  //printf("sd.scaler->mySd.jout=%d\n",scaler.mySd.jout);
  
  // OPTIONS PARSING
  // BMA TODO: max_trials, lower/upper bounds (use toVec helper),
  //           correlation_lengths, optimizationMethod
  ParamMap::const_iterator param_it;


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
  iAnchorPoint=0;
  ifHaveAnchorPoint=false;
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
    cerr << "Your options are to\n(A) specify both the upper and lower, or\n(B) specify neither the upper nor lower,\nbounds of the domain of the Kriging Model\n";
    assert(if_user_specified_lower_bounds==if_user_specified_upper_bounds);
  }
  
  if(if_user_specified_lower_bounds==true) {
    for(int ivarsr=0; ivarsr<numVarsr; ++ivarsr) 
      if(!(min_max_xr(0,ivarsr)<=min_max_xr(1,ivarsr))) {
	cerr << "The lower bound of the domain of the Kriging Model must be less than or equal to the upper bound of the domain of the Kriging Model\n";
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
  
  //printf("KrigingModel constructor should have just written out domain bounds\n");
  
  param_it = params.find("dimension_groups");
  if (param_it != params.end() && param_it->second.size() > 0) {
    MtxInt dim_groups(1,numVarsr);
    assert(!(dim_groups.putRows(param_it->second,0)));
    sdBuild.setDimGroups(dim_groups);
  }
  
  scaler.scaleToDefault(); //scale outputs to -0.5<=Y<=0.5 and scale real inputs to volume 1 hyper-rectangle centered at 0 if real iputs dimensions are locked or the unit hypercube centered at 0 if no dimensions are locked.  The scaling is done to let us define the feasible region simply (done in create);

  sdBuild.getY(Yall);  

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
  else if(optimizationMethod.compare("global")==0)
    maxTrials = 10000;
  else if(optimizationMethod.compare("global_local")==0) {
    maxTrials = 10000; //ensure it has non-zero as a fail safe but this shouldn't be used
    maxTrialsGlobal = 500;
    maxTrialsLocal = 20;
  }
  else{ //error checking the input
    cerr << "KrigingModel() unknown optimization_method [" << optimizationMethod << "]  aborting\n";
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
    cerr << "Local optimization is the only optimization method for Kriging that uses the \"num_starts\" key word. Check your input file for errors.\n";
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
      cerr << "You can't both \n (A) use the global optimization method to choose, and \n (B) directly specify \n correlation lengths for the Kriging model.\n";
      assert(optimizationMethod.compare("global")!=0);
    }
    else if(optimizationMethod.compare("sampling")==0) {
      // this is only the default number of samples/maxTrials, the user can 
      // still overide this below
      maxTrials+=1;
    }
    
    natLogCorrLen.newSize(1,numVarsr); //allocate space 
    
    //read the correlation lengths in from the string
    assert(!(natLogCorrLen.putRows(param_it->second,0)));
    // "natLogCorrLen" currently holds the unscaled correlation LENGTHS, not 
    // the natural log of the scaled correlation length, we need to fix that
    // but first we need to check the input for errors
    for(int ivarsr=0; ivarsr<numVarsr; ++ivarsr) 
      if(!(natLogCorrLen(0,ivarsr)>0.0)) {
	cerr << "For the Kriging Model, correlation lengths must be strictly positive\n.";
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

  //cout << "order=" << polyOrder << "\n";

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
  numConFunc=1;
  
  //convert to the Dakota bitflag convention for derivative orders
  int num_analytic_obj_ders_in=0;
  int num_analytic_con_ders_in=0;
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

  //nug=(2*nTrend+1.0)/maxCondNum;
  //nug=numPoints/(maxCondNum-1.0); //this is guaranteed to make R non-singular it's a little overkill though

  ifChooseNug = false;
  param_it = params.find("find_nugget");
  if (param_it != params.end() && param_it->second.size() > 0)
    ifChooseNug = true; 

  //ifChooseNug = true ; std::cout << "ifChooseNug=" << ifChooseNug << "\n";

  // fixed value for now
  maxChooseNug = 0.1;
  nug = 0.0; //default

  nuggetFormula=0;
  param_it = params.find("nugget_formula");
  if (param_it != params.end() && param_it->second.size() > 0) {
    if(ifChooseNug==true) {
      cerr << "You can't both auto-select a nugget and use a preset formula" << endl;
      assert(ifChooseNug==false);
    }
    nuggetFormula=std::atoi(param_it->second.c_str()); 
    if(nuggetFormula!=0) {
      switch(nuggetFormula) {
      case 1:
	nug=(2*numTrend(polyOrderRequested,0)+1.0)/maxCondNum;
	break;
      case 2:
	nug=2*numPoints/maxCondNum;
	break;
      default:
	cerr << "nugget_formula =" << nuggetFormula << " is not one of the available preset nugget formulas." << endl;
	assert(false);
      }
    }
  }

  param_it = params.find("nugget");
  if (param_it != params.end() && param_it->second.size() > 0) {
    if(!((nuggetFormula==0)&&(ifChooseNug==false))) {
      cerr << "You can do at most 1 of the following (A) auto-select the nugget (approximately the minimum needed to satisfy the condition number bound) (B) use one of the preset nugget formulas (C) directly specify a nugget.  The default is not to use a nugget at all (i.e. use a nugget of zero)." << endl;
      assert((nuggetFormula==0)&&(ifChooseNug==false));
    }
    nug = std::atof(param_it->second.c_str()); 
    if(!(nug >= 0.0)) {
      cerr << "The nugget must be greater than or equal to zero." << endl;
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
  //lrp.optimize_by_picking_best_guess(EulAng); 
  //EulAng.reshape(nchoosek(numVarsr, 2),1);
  EulAng.newSize(nchoosek(numVarsr, 2),1); 
  EulAng.zero();
  gen_rot_mat(Rot, EulAng, numVarsr);
  eval_trend_fn(Gall, Poly, Rot, XR);
  //LinearRegressionModel::evalBasis(G,poly,Rot,XR);

  gen_Z_matrix();  

  //printf("completed the right KrigingModel constructor\n", stdout); fflush(stdout);
}

void KrigingModel::create()
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

  //  if(numPoints<=2*numVarsr)
  //  aveDistBetweenPts=1.0;
  //else
  //  aveDistBetweenPts=std::pow(numPoints-2*numVarsr,-1.0/numVarsr);

  //printf("aveDistBetweenPts=%12.6g\n",aveDistBetweenPts);
  /* Gaussian Process error model has about ~5% confidence (2 std devs away) 
     in what points 4 neighbors away have to say. If points are correlated 
     well at even greater distances then either 
     * that same information will be contained in nearby points OR
     * you shouldn't be using a Gaussian process error model
     KRD */
  double max_corr_length = aveDistBetweenPts*8.0; //*2.0; 

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

  //printf("got to yada yada\n"); fflush(stdout);
  ///set the bounds and the initial iterates
  if(ifUserSpecifiedCorrLengths==true) {
    //printf("says that the user specified correlation lengths\n"); fflush(stdout);
    // the first guess is what the user told us he/she wanted to use
    for(int jvar=0; jvar<numVarsr; ++jvar) {
      opt.lower_bound(jvar, minNatLogCorrLen);
      opt.upper_bound(jvar, maxNatLogCorrLen);
      //double temp_double=natLogCorrLen(0,jvar);
      //printf("KrigingModel::create() jvar=%d temp_double=%g\n",jvar,temp_double);
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
  //bins are the endpoints of a randomly rotated axis
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
      cerr << "KrigingModel:create() unknown optimization_method [" << optimizationMethod << "]  aborting\n";
      assert(false);
    }
    natLogCorrLen = opt.best_point();
  }

  correlations.newSize(1,numVarsr);
  //MtxDbl *mtxdblptr=&natLogCorrLen;
  //int nROWS=natLogCorrLen.getNRows();
  //int nCOLS=natLogCorrLen.getNCols();
  //printf("mtxdblptr=%d nROWS=%d nCOLS=%d",mtxdblptr,nROWS,nCOLS);
  //fflush(stdout);
  //printf("\n");
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

  temp.clear();
  temp2.clear();

  //variables whose values needed to be retained between sequential call to masterObjectiveAndConstraints for precompute and store strategy to work
  prevObjDerMode=prevConDerMode=0;
  prevTheta.clear(); //row vector 
  Z.clear(); //matrix
  Ztheta.clear(); //vector
  R.clear(); //matrix
  G.clear(); //matrix
  Yall.clear();
  Gall.clear();
  rcondDblWork.clear();
  rcondIntWork.clear();
  sumAbsColR.clear();
  lapackRcondR.clear();
  oneNormR.clear();

  con.clear(); //vector
  gradObj.clear(); //vector
  gradCon.clear(); //matrix
  hessObj.clear(); //matrix used to compute hessObj
}

std::string KrigingModel::model_summary_string() const {
  MtxDbl temp_out_corr_lengths(1,numVarsr);
  for(int i=0; i<numVarsr; ++i) 
    temp_out_corr_lengths(0,i)=std::sqrt(0.5/correlations(0,i));
  scaler.unScaleXrDist(temp_out_corr_lengths);
  
  //printf("numPoints=%d numTrend=%d numEqnKeep=%d\n",numPoints,numTrend(polyOrder,0),numEqnKeep);


  std::ostringstream oss;
  oss << "--- Surfpack Kriging Diagnostics ---\n";
  oss << "KM: #pts=" << numPoints << "; used " << numEqnKeep << "/" << numEqnAvail << " pts; Correlation lengths=(" << temp_out_corr_lengths(0,0);
  for(int i=1; i<numVarsr; ++i)
    oss << ", " << temp_out_corr_lengths(0,i);
  oss << "); unadjusted variance=" << estVarianceMLE * scaler.unScaleFactorVarY() << "; log(likelihood)=" << likelihood << "; the trend is a ";
  if(polyOrder>1) {
    if(ifReducedPoly==true)
      oss << "reduced_";
    else oss <<"full ";
  }
  oss << "polynomial of order=" << polyOrder << "(order " << polyOrderRequested << 
    " was desired); rcond(R)=" << rcondR << "; rcond(Gtran_Rinv_G)=" << rcond_Gtran_Rinv_G 
      << "; nugget=" << nug << ".\n";
  //printf("mod_sum_str nug=%22.16g\n",nug);

  oss << "Beta= (" << betaHat(0,0);
  for(int i=1; i<nTrend; ++i)
    oss << "," << betaHat(i,0);
  oss << ")\n";
  oss << "------------------------------------\n";
  return (oss.str());  
}

// BMA TODO: combine these two functions?

/// evaluate (y) the Kriging Model at a single point (xr)
double KrigingModel::evaluate(const MtxDbl& xr) const
{
  double singular_y;
  if(scaler.isYSingular(0,singular_y)) {
    //you wouldn't want to do this for gradient based Kriging
    //if gradients of y were used as inputs
    return singular_y;
  }

  //assert( (numVarsr == xr.getNCols()) && (xr.getNRows() == 1) );
  MtxDbl g(1, nTrend), r(1, numRowsR);

  /*
  printf("double evaluate()\n");
  printf("xr=[%20.14g", xr(0,0));
  for(int i=1; i<numVarsr; ++i)
    printf(", %20.14g",xr(0,i));
    printf("]\n");
  */
  
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
      
  //if(nug>0.0)
  //apply_nugget_eval(r);
  
  double y = dot_product(g, betaHat) + dot_product(r, rhs);

  //double yus=scaler.unScaleYOther(y);
  //printf("] y=%g\n",yus);
  //printf("y=%g yunscaled=%g\n",y,yus);
  //return yus;

  return (scaler.unScaleYOther(y));
}


/// evaluate (y) the Kriging Model at a collection of points (xr)
MtxDbl& KrigingModel::evaluate(MtxDbl& y, const MtxDbl& xr) const
{
  int nrowsxr = xr.getNRows();
  //printf("nrowsxr=%d nvarsrxr=%d",nrowsxr,xr.getNCols());

  y.newSize(nrowsxr, 1);
  double singular_y;
  if(scaler.isYSingular(0,singular_y)) {
    //you wouldn't want to do this for gradient based Kriging
    //if gradients of y were used to build model
    for(int i=0; i<nrowsxr; ++i)
      y(i,0)=singular_y;
    return y;
  }
  
  //assert( (numVarsr == xr.getNCols()) && (xr.getNRows() == 1) );
  MtxDbl g(nrowsxr, nTrend), r(nrowsxr, numRowsR);

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

  //  if(nug>0.0)
  //apply_nugget_eval(r);
  
  //y=0.0*y+1.0*g*betaHat => y=g*betaHat
  matrix_mult(y, g , betaHat, 0.0, 1.0);
  
  //y=1.0*y+1.0*r*rhs where rhs=R^-1*(Y-G(XR)*betaHat), initial y=g*betaHat => y=g*betaHat+r*rhs
  matrix_mult(y, r, rhs    , 1.0, 1.0);
  
  scaler.unScaleYOther(y);
  //printf("y is correct for ValidateMain because it isn't being unscaled\n");

  return y;
}

MtxDbl& KrigingModel::evaluate_d1y(MtxDbl& d1y, const MtxDbl& xr) const
{
  int nrowsxr = xr.getNRows();
  d1y.newSize(nrowsxr, numVarsr);
  double singular_y;
  if(scaler.isYSingular(0,singular_y)) {
    //you wouldn't want to do this for gradient based Kriging
    //if gradients of y were used to build model
    d1y.zero();
    return d1y;
  }

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
  //apply_nugget_eval(r);
  MtxDbl d1r(nrowsxr,numRowsR);
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

    dcorrelation_matrix_dxI(d1r, r, xr_scaled, ivar);
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

MtxDbl& KrigingModel::evaluate_d2y(MtxDbl& d2y, const MtxDbl& xr) const
{
  int nrowsxr=xr.getNRows();
  int nder=num_multi_dim_poly_coef(numVarsr,-2);
  d2y.newSize(nrowsxr,nder);
  double singular_y;
  if(scaler.isYSingular(0,singular_y)) {
    //you wouldn't want to do this for gradient based Kriging
    //if gradients of y were used as inputs
    d2y.zero();
    return d2y;
  }

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
  //apply_nugget_eval(r);
  MtxDbl d1r(nrowsxr,numRowsR);
  MtxDbl d2r(nrowsxr,numRowsR);
  MtxDbl temp_vec(nrowsxr,1);

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
      dcorrelation_matrix_dxI(d1r, r, xr_scaled, ivar);
    }

    //find the second dimension we are taking a first derivative of
    if(der(ider,ivar)==2)
      jvar=ivar;
    else
      for(jvar=ivar+1; jvar<numVarsr; ++jvar)
	if(der(ider,jvar)>0)
	  break;
    
    //dcorrelation_matrix_dxI(d2r, d1r, xr_scaled, jvar);
    d2correlation_matrix_dxIdxK(d2r, d1r, r, xr_scaled, ivar, jvar);
    
    matrix_mult(temp_vec,d2r,rhs);

    for(int ipt=0; ipt<nrowsxr; ++ipt)
      d2y(ipt,ider)=(d2y(ipt,ider)+temp_vec(ipt,0))*d2y_unscale_factor;
  }

  return d2y;
}



/// matrix Ops evaluation of adjusted variance at a single point
double KrigingModel::eval_variance(const MtxDbl& xr) const
{
  double singular_y;
  if(scaler.isYSingular(0,singular_y)) {
    //you wouldn't want to do this for gradient based Kriging
    //if gradients of y were used as inputs
    //printf("NKM eval_variance: y is singular\n");
    return 0;
  }

  //assert( (numVarsr == xr.getNCols()) && (xr.getNRows() == 1) );
  MtxDbl g_minus_r_Rinv_G(1, nTrend), r(1, numRowsR);

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

  MtxDbl tempa(numRowsR,1);
  MtxDbl tempb(nTrend,1);

  //if(nug>0.0) 
  //apply_nugget_eval(r);

  solve_after_Chol_fact(tempa,RChol,r,'T');
  
  matrix_mult(g_minus_r_Rinv_G,r,Rinv_G,1.0,-1.0);
  solve_after_Chol_fact(tempb,Gtran_Rinv_G_Chol,g_minus_r_Rinv_G,'T');

  double unscale_factor_vary=scaler.unScaleFactorVarY();
  double adj_var=estVarianceMLE*unscale_factor_vary*
    (1.0-dot_product(tempa,r)+dot_product(tempb,g_minus_r_Rinv_G));
  //if(!(adj_var>0.0)) {
  //printf("adj_var=%g estVarianceMLE=%g rcondR=%g unscale_factor_vary=%g\n",adj_var,estVarianceMLE,rcondR,unscale_factor_vary); 
  //fflush(stdout);
  //}
  adj_var=std::fabs(adj_var); //hack to handle negative zero variance (machine precision round off error)
  if(adj_var<0.0) {
    printf("NKM setting adj_var to zero adj_var=%g unadj_var=%g rcondR=%g\n",adj_var,estVarianceMLE*unscale_factor_vary,rcondR); 
    adj_var=0.0;
  }
  else if(adj_var==0.0)
    printf("NKM adj_var is zero =%g\n",adj_var);
  else if(!(adj_var>=0.0))
    printf("double NKM_KrigingModel::eval_variance(...) adj_var=nan rcondR=%g\n",rcondR);

  return adj_var;
}

/// matrix Ops (as much as possible with BLAS and LAPACK) evaluation of adjusted variance for a collection of points... The MATLAB would be estVarianceMLE*(1-sum((r/R).*r,2)+sum((g_minus_r_Rinv_G/(Gtran_Rinv_G)).*g_minus_r_Rinv_G,2) unfortunately there's not a convenient way to do it with BLAS & LAPACK
MtxDbl& KrigingModel:: eval_variance(MtxDbl& adj_var, const MtxDbl& xr) const
{
  int nrowsxr=xr.getNRows();
  adj_var.newSize(nrowsxr,1);

  double singular_y;
  if(scaler.isYSingular(0,singular_y)) {
    //you wouldn't want to do this for gradient based Kriging
    //if gradients of y were used as inputs
    adj_var.zero();
    return adj_var;
  }

  //assert( (numVarsr == xr.getNCols()) && (xr.getNRows() == 1) );
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

  //if(nug>0.0)
  //apply_nugget_eval(r);
  
  solve_after_Chol_fact(tempa,RChol,r,'T');

  matrix_mult(g_minus_r_Rinv_G,r,Rinv_G,1.0,-1.0);
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
      printf("MtxDbl& NKM_KrigingModel::eval_variance(...) adj_var(%d)=nan rcondR=%g\n",i,rcondR);
	

  return adj_var;
}


/*
VecDbl KrigingModel::gradient(const VecDbl& x) const
{
  assert(!x.empty());
  assert(x.size()+1==betaHat.size()); //true for linear trend function; KRD added
  cout << "IN gradient x[0] = " << x[0] << endl;
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

std::string KrigingModel::asString() const
{
  std::ostringstream os;
  os << "\ncenters:\n" << bs.asString() << "rhs: ";
  copy(rhs.begin(),rhs.end(),std::ostream_iterator<double>(os," "));
  os << "\n";
  return os.str();
}
*/


#ifdef __DONT_USE__
/* this function scales the input matrix, r, in place by 1.0/(1.0+nug), 
    where nug is the nugget (modifying the correlation matrix by the 
    inclusion of a nugget causes the KrigingModel to smooth data, i.e. 
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

MtxDbl& KrigingModel::apply_nugget_eval(MtxDbl& r) const {
  if(!(nug>0.0))
    return r;

  //printf("apply_nugget_eval\n");

  //int nelem=r.getNElems();
  int nrowsr=r.getNRows();
  double temp_dbl=1.0/(1.0+nug);

  for(int j=0; j<numRowsR; ++j)
    for(int i=0; i<nrowsr; ++i)
      r(i,j)*=temp_dbl;
  //for(int ij=0; ij<nelem; ++ij) 
  //r(ij)*=temp_dbl;

  //printf("apply_nugget_eval temp=%g\n",temp);

  return r;
}
#endif


/** set R=(1.0+nug)^-1*(R+nug*I), where R is the correlation matrix for the
    data that the model is built from.  This function preserved the original 
    R.  Modifying the correlation matrix by the inclusion of a nugget causes 
    the KrigingModel to smooth the data, i.e. approximate it rather than 
    interpolate it, it can also be used to fix an ill conditioned correlation 
    matrix.  The convention is that capital matrices are for data the model 
    is built from, lower case matrices are for arbitrary points to evaluate 
    the model at */
void KrigingModel::apply_nugget_build() {
  if(!(nug>0.0)) return;
  //printf("applying nugget=%22.16g\n",nug);

  int nrowsR=R.getNRows();
  //assert(nrowsR==R.getNCols());

  //if I do it this way I don't need to call apply_nugget_eval()
  for(int i=0; i<nrowsR; ++i)
    R(i,i)+=nug;

  /*

  double temp_dbl=1.0/(1.0+nug);
  for(int j=0; j<nrowsR; ++j)
    for(int i=0; i<nrowsR; ++i)
      R(i,j)*=temp_dbl;

  //switching to matrices that can have a different (smaller) apparent than actual (memory footprint size) and uses index math rather than a  array of pointers to column leaders (which was committed to the repository yesteday, 2011.05.31), means that division and mod operators are needed for single index lookup, which makes single index look up slower than 2 index look up, so switching to 2 index lookup on 2011.06.01)
  //int nelemsR=nrowsR*nrowsR;
  //int ij;
  //for(ij=0; ij<nelemsR; ++ij)
  //R(ij)*=temp_dbl;
  
  //the "paranoid" part of my mind wonders if there would be less round off
  //error if I added the nugget to the diagonal BEFORE scaling, the pragmatic
  //part of my mind says it shouldn't matter and doing it this way is faster
  temp_dbl*=nug;
  for(int i=0; i<nrowsR; ++i)
    R(i,i)+=temp_dbl; 
  */
  return;
}

/** r (lower case r) is the correlation matrix between the
    interpolation points and data points, it used to EVALUATE but not
    construct the emulator's Gaussian process error model
    i.e. E(y(xr)|Y(XR))=g(xr)*betaHat+r*R^-1*eps where
    eps=(Y-G(XR)*betaHat), to be more specific
    r(i,j)=r(xr(i,:),XR(j,:))=exp(sum_k -theta(0,k)*(xr(i,k)-XR(j,k))^2) KRD
    wrote this */
MtxDbl& KrigingModel::correlation_matrix(MtxDbl& r, const MtxDbl& xr) const
{
  //int nrowsXR=XR.getNRows(); //data points used to build model
  int nrowsxr=xr.getNRows(); //points at which we are evalutating the model
  //assert(xr.getNCols()==numVarsr);
  r.newSize(nrowsxr,numRowsR);
  int i,j,k, jeqn;
  double temp_double;

  /*
  printf("**********************************************************\n");
  printf("**********************************************************\n");
  for(i=0; i<nrowsXR; ++i) {
    printf("XR(%d,:)={",i);
    for(j=0; j<numVarsr; ++j) printf(" %g",XR(i,j));
    printf(" }; Y(%d)=%g\n",i,Y(i,0));
  }
  printf("**********************************************************\n");
  printf("**********************************************************\n");
  */
  
  if(numVarsr==1) {
    //optimized for when there is only 1 output variable
    for(j=0; j<numRowsR; ++j) {
      jeqn=iEqnKeep(j,0);
      for(i=0; i<nrowsxr; ++i) {
	temp_double=xr(i,0)-XR(jeqn,0);
	r(i,j)=std::exp(-correlations(0,0)*temp_double*temp_double); 
      }
    }
  } else if(nrowsxr==1) {
    //"optimized" for when there is only 1 evaluation point (Save loops, save dereferences)
    
    //k=0 was pulled out from below to avoid doing an extra loop just to 
    //initialize all of r to zero
    for(j=0; j<numRowsR; ++j) {
      temp_double=xr(0,0)-XR(iEqnKeep(j,0),0);
      r(0,j)=-correlations(0,0)*temp_double*temp_double; //=- is correct
    }
  
    //all but first and last k
    for(k=1;k<numVarsr-1;++k)
      for(j=0; j<numRowsR; ++j) {
	temp_double=xr(0,k)-XR(iEqnKeep(j,0),k);
	r(0,j)-=correlations(0,k)*temp_double*temp_double; //-= is correct
      }
  
    //this value of k was pulled out of above to save doing an extra loop 
    //for just the exp() operation
    k=numVarsr-1; 
    for(j=0; j<numRowsR; ++j) {
      temp_double=xr(0,k)-XR(iEqnKeep(j,0),k);
      r(0,j)=std::exp(r(0,j)-correlations(0,k)*temp_double*temp_double); 
    }
  } else {
    //"optimized" for when there is more than 1 output variable
    
    //k=0 was pulled out from below to avoid doing an extra loop just to 
    //initialize all of r to zero
    for(j=0; j<numRowsR; ++j) {
      jeqn=iEqnKeep(j,0);
      for(i=0; i<nrowsxr; ++i) {
	temp_double=xr(i,0)-XR(jeqn,0);
	r(i,j)=-correlations(0,0)*temp_double*temp_double; //=- is correct
      }
    }
    //all but first and last k
    for(k=1;k<numVarsr-1;++k)
      for(j=0; j<numRowsR; ++j) {
	jeqn=iEqnKeep(j,0);
	for(i=0; i<nrowsxr; ++i) {
	  temp_double=xr(i,k)-XR(jeqn,k);
	  r(i,j)-=correlations(0,k)*temp_double*temp_double; //-= is correct
	}
      }

    //this value of k was pulled out of above to save doing an extra loop 
    //for just the exp() operation
    k=numVarsr-1; 
    for(j=0; j<numRowsR; ++j) {
      jeqn=iEqnKeep(j,0);
      for(i=0; i<nrowsxr; ++i) {
	temp_double=xr(i,k)-XR(jeqn,k);
	r(i,j)=std::exp(r(i,j)-correlations(0,k)*temp_double*temp_double); 
      }
    }
  }
  return r;
}

///k is the variable/dimension not the point
MtxDbl& KrigingModel::dcorrelation_matrix_dxI(MtxDbl& dr, const MtxDbl& r, 
					      const MtxDbl& xr, int Ider) const
{
  int nrowsxr=xr.getNRows();
  assert((r.getNRows()==nrowsxr)&&(r.getNCols()==numRowsR)&&
	 (xr.getNCols()==numVarsr)&&(0<=Ider)&&(Ider<numVarsr));
  dr.newSize(nrowsxr,numRowsR);

  int jeqn;
  double temp_dbl=-2.0*correlations(0,Ider);
  for(int j=0; j<numRowsR; ++j) {
    jeqn=iEqnKeep(j,0);
    for(int ipt=0; ipt<nrowsxr; ++ipt)
      dr(ipt,j)=temp_dbl*r(ipt,j)*(xr(ipt,Ider)-XR(jeqn,Ider));
  }
  return dr;
}

MtxDbl& KrigingModel::d2correlation_matrix_dxIdxK(MtxDbl& d2r, const MtxDbl& drI, const MtxDbl& r, const MtxDbl& xr, int Ider, int Kder) const
{
  //int nrowsXR=XR.getNRows(); //data points used to build model
  int nrowsxr=xr.getNRows(); //points at which we are evalutating the model
  d2r.newSize(nrowsxr,numRowsR);
///k is the variable/dimension not the point

  assert((r.getNRows()==nrowsxr)&&(r.getNCols()==numRowsR)&&
	 (xr.getNCols()==numVarsr)&&(0<=Kder)&&(Kder<numVarsr));

  double neg_two_theta_K=-2.0*correlations(0,Kder);
  int jeqn;
  if(Ider==Kder) 
    for(int j=0; j<numRowsR; ++j) {
      jeqn=iEqnKeep(j,0);
      for(int ipt=0; ipt<nrowsxr; ++ipt)
	d2r(ipt,j)=neg_two_theta_K*((xr(ipt,Kder)-XR(jeqn,Kder))*drI(ipt,j)+r(ipt,j));
    }
  else
    for(int j=0; j<numRowsR; ++j) {
      jeqn=iEqnKeep(j,0);
      for(int ipt=0; ipt<nrowsxr; ++ipt)
	d2r(ipt,j)=neg_two_theta_K*(xr(ipt,Kder)-XR(jeqn,Kder))*drI(ipt,j);
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
    is produced by KrigingModel::gen_Z_matrix()     KRD wrote this */
void KrigingModel::correlation_matrix(const MtxDbl& theta)
{
  int nrowsZ=Z.getNRows();
  //printf("nrowsZ=%d; numPoints=%d; ''half'' numPoints^2=%d; numVarsr=%d; theta.getNCols()=%d\n",
  //	 nrowsZ,numPoints,nchoosek(numPoints,2),numVarsr,theta.getNCols());
  //fflush(stdout);
  //assert((nrowsZ==nchoosek(numPoints,2))&&(numVarsr==theta.getNCols()));

  Ztheta.newSize(nrowsZ,1);
  matrix_mult(Ztheta,Z,theta,0.0,1.0,'N','T');
  numRowsR=numEqnAvail; //numEqnAvail=numPoints for regular (not derivative enhanced Kriging)
  R.newSize(numPoints,numPoints);

  double Rij_temp;
  int ij=0;
  for(int j=0; j<numPoints-1; ++j) {
    R(j,j)=1.0;
    for(int i=j+1; i<numPoints; ++i, ++ij) {
      Rij_temp=std::exp(Ztheta(ij,0));
      R(i,j)=Rij_temp;
      R(j,i)=Rij_temp;
    }
  }
  R(numPoints-1,numPoints-1)=1.0;

  //  FILE *fp=fopen("km_Rmat_check.txt","w");
  //  for(int i=0; i<numPoints; ++i) {
  //    fprintf(fp,"%-12.6g", R(i,0));
  //    for(int j=1; j<numPoints; ++j) 
  //      fprintf(fp," %-12.6g", R(i,j));     
  //    fprintf(fp,"\n");
  //  }
  //  fclose(fp);

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
MtxDbl& KrigingModel::gen_Z_matrix()
{
  int nrowsXR=XR.getNRows();
  int ncolsXR=XR.getNCols();
  int nrowsZ=nchoosek(nrowsXR,2);

  Z.newSize(nrowsZ,ncolsXR);

  register double mult_term;
  //register int ijk=0;
  //double *Z_ptr=Z.ptr(0); //done for speed, undone for robustness because matrices can now have fewer apparent rows than actual rows (it shouldn't happen for Kriging that isn't gradient enhanced, but future mods might change that) 2011.06.01
  const double *XR_k_ptr; //done for speed
  for(int k=0; k<ncolsXR; k++) {
    XR_k_ptr=XR.ptr(0,k);
    int ij=0;
    for(int j=0; j<nrowsXR-1; ++j)
      for(int i=j+1; i<nrowsXR; ++i, ++ij) {
	mult_term=XR_k_ptr[i]-XR_k_ptr[j];
	//Z_ptr[ijk++]=-mult_term*mult_term;
	Z(ij,k)=-mult_term*mult_term;
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
void KrigingModel::masterObjectiveAndConstraints(const MtxDbl& theta, int obj_der_mode, int con_der_mode)
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
  assert((maxObjDerMode<=7)&&(maxConDerMode<=3)&&
	 (0<=obj_der_mode)&&(obj_der_mode<=maxObjDerMode)&&
	 (0<=con_der_mode)&&(con_der_mode<=maxConDerMode)&&
	 ((1<=obj_der_mode)||(1<=con_der_mode)));

  //next 5 or 6 must be copied out (can't allow external code to change these or have no guarantee that we can reuse "previous" values)
  //member variable: double obj;
  //member variable: MtxDbl con(numConFunc);

  //these 2 are set by choice of optimizer and are not changed inside this function
  //member variable: int maxObjDerMode;
  //member variable: int maxConDetMode;

  //these are private and their values need to be retained between sequential calls to masterObjectiveAndConstraints, that is no other function (other than the create) can access them, these get deallocated at the end of create
  //member variable: int prevObjDerMode
  //member variable: int prevConDerMode
  //member variable: MtxDbl prevTheta(1,numTheta)
  //member variable: MtxDbl Z(numPoints*numPoints,numTheta)
  //member variable: MtxDbl R(numPoints,numPoints)
  //member variable: MtxDbl G(numPoints,nTrend)

  //keep these around after emulator creation so we can evaluate the emulator
  //member variable: MtxDbl rhs(numPoints)
  //member variable: MtxDbl RChol(numPoints,numPoints)
  //member variable: MtxDbl Rinv(numPoints,numPoints) //needed to eval integral of adjusted variance, capability not yet added
  //member variable: MtxDbl Rinv_G(numPoints,nTrend)
  //member variable: MtxDbl Gtran_Rinv_G_Chol(nTrend,nTrend)

  //if theta was the same as the last time we called this function than we can reuse some of the things we calculated last time
  
  int i;

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
  if((prevObjDerMode==0)&&(prevConDerMode==0))
    for(i=0; i<numTheta; ++i) 
      prevTheta(0,i)=theta(0,i); 

  int chol_info;

  if((prevObjDerMode==0)&&(prevConDerMode==0)) {
    R.newSize(numEqnAvail,numEqnAvail);
    correlation_matrix(theta); //fills member variable R as exp(Z*theta) where Z is a member variable
    //apply_nugget_build(); //modify R by nug in place
  }

  if((prevObjDerMode==0)&&((1<=obj_der_mode)||(1<=con_der_mode))) {


    //perform LU decomposition of R and calculate the determinant of R, replaced LU with Cholesky
    //http://en.wikipedia.org/wiki/Determinant#Determinant_from_LU_decomposition

    if(ifChooseNug==true)
      nuggetSelectingCholR();
    else
      equationSelectingCholR();
    double min_allowed_rcond=1.0/maxCondNum;
    if((rcondR<=min_allowed_rcond)||(numRowsR<=numTrend(polyOrder,0))) {
      //nTrend=numTrend(polyOrder,0);
      printf("singular correlation matrix rcondR=%g numRowsR=%d numTrend=%d\n",
	     rcondR,numRowsR,numTrend(polyOrder,0));
      obj=HUGE_VAL;
      con.newSize(numConFunc,1);
      for(int i=0; i<numConFunc; ++i)
	con(i,0)=1.0; //say the constraints are violated
      return;
    }

    //RChol.copy(R);
    //chol_info=0;
    //Chol_fact(RChol,chol_info,rcondR); //preconditioned Cholesky, when Kriging is not gradient enhaced R won't need preconditiong since it has all 1's on the diagonals
    //printf("Chol\n");
    //assert(chol_info==0);  //here for debug, decide what to do about it later

    double log_determinant_R = 0.0; //need to do this to avoid underflow error for large numbers of points, log(0)=-inf
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
    solve_after_Chol_fact(Rinv_G,RChol,G);
    //solve_after_LU_fact(Rinv_G,RLU,ipvt_RLU,G,'N','N'); //O(N^3) ops

    Gtran_Rinv_G_Chol.newSize(nTrend,nTrend);
    matrix_mult(Gtran_Rinv_G_Chol,G,Rinv_G,0.0,1.0,'T','N');
    //double rcond_Gtran_Rinv_G;
    //Chol_fact(Gtran_Rinv_G_Chol,chol_info,rcond_Gtran_Rinv_G);
    Chol_fact_workspace(Gtran_Rinv_G_Chol,Gtran_Rinv_G_Chol_Scale,Gtran_Rinv_G_Chol_DblWork,Gtran_Rinv_G_Chol_IntWork,chol_info,rcond_Gtran_Rinv_G);

    double log_determinant_Gtran_Rinv_G=0.0;
    for (int itrend = 0; itrend < nTrend; ++itrend)
      log_determinant_Gtran_Rinv_G += std::log(Gtran_Rinv_G_Chol(itrend,itrend)); 
    log_determinant_Gtran_Rinv_G *= 2.0; //only for Cholesky factorization of R 
    //if(~(chol_info==0)) assert(chol_info==0);  //for debug, do something else for production

    temp.newSize(nTrend,1);
    matrix_mult(temp, Rinv_G, Y, 0.0, 1.0, 'T', 'N');
    betaHat.newSize(nTrend,1);


    solve_after_Chol_fact(betaHat,Gtran_Rinv_G_Chol,temp); //O(nTrend^2) ops
    temp2.copy(Y); //this will be eps=epsilon=Y-G(XR)*betaHat, but use 
    //variable temp2 because we would only need the variable "eps" for these 
    //5 lines of code (not counting comments) and we want to save space, 
    //afterwards we will only need R^-1*eps which is stored in "rhs"
    matrix_mult(temp2, G, betaHat, 1.0, -1.0, 'N', 'N'); //eps=Y-G(XR)*betaHat
    rhs.newSize(numRowsR,1);
    solve_after_Chol_fact(rhs,RChol,temp2);




    //it's actually the log likelihood, which we want to maximize
    //likelihood = -0.5*(numPoints*(std::log(4.0*std::acos(0.0))+std::log(estVarianceMLE)+1)
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

    //the "per point" log(likelihood)
    likelihood = -0.5*(std::log(estVarianceMLE)+log_determinant_R/numRowsR); 
#endif

    //if(likelihood>=DBL_MAX)
    //printf("[estVarianceMLE=%g determinant_R=%g]",estVarianceMLE,determinant_R);

    //the objective function being MINIMIZED is the negative of the log 
    //likelihood (on a per point basis so numbers will be comparable 
    //regardless of how many points there are)
    obj=-likelihood;  
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
    else{
      printf("Error:  the only current option for constraint type is \"rcond\"\n");
      assert(constraintType.compare("rcond")==0);
    }

    prevConDerMode=1; //increase prevConDerMode to current value
    if((con_der_mode==1)&&(obj_der_mode<=prevObjDerMode)) {
      //we have everything we need so exit early
      return;
    }
  }
  
  return;
}


void KrigingModel::getRandGuess(MtxDbl& guess) const
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
    //guess(0,j) = 1.0/(2.0*corr_length*corr_length); 
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
MtxDbl& KrigingModel::makeGuessFeasible(MtxDbl& nat_log_corr_len, OptimizationProblem *opt) {
  int k;
  MtxDbl theta(1,numTheta);
  for(k=0; k<numTheta; ++k)
    theta(0,k)=0.5*std::exp(-2.0*nat_log_corr_len(0,k));

  R.newSize(numEqnAvail,numEqnAvail);
  
  correlation_matrix(theta); //assigns to member variable R

  //if((ifChooseNug==true)||(nug<0.0))
  //nug=0.0;

  //need to revisit this to allowed specified nug
  //apply_nugget_build();
  //RChol.copy(R);
  double best_rcond;
    if(ifChooseNug==true)
      nuggetSelectingCholR();
    else
      equationSelectingCholR();
    //equationSelectingCholR();
  if(numRowsR<=numTrend(polyOrder,0)) 
    best_rcond=0;
  else
    best_rcond=rcondR;
  //Chol_fact(RChol,chol_info,best_rcond);
  //if(~(chol_info==0)) assert(chol_info==0); //for debug, do something different for production
  //RLU.copy(R);
  //LU_fact(RLU,ipvt_RLU);
  //double best_rcond=rcond_after_LU_fact(R,RLU);
  //double best_rcond=rcond_after_Chol_fact(R,RChol);
  double rcond;
  MtxDbl guess(1,numTheta);
  MtxDbl guess_theta(1,numTheta);
  int iguess=0;
  while((1.0-best_rcond*maxCondNum>0.0)&&(iguess<50)) {
    //while((1.0-best_rcond*3.0*maxCondNum<0.0)&&(iguess<50)) {
    iguess++;
    getRandGuess(guess);

    //convert guess from nat_log_corr_len to theta
    for(k=0; k<numTheta; ++k)
      guess_theta(0,k)=0.5*std::exp(-2.0*guess(0,k));

    correlation_matrix(guess_theta); //assigns to member variable R

    //need to revisit this to allowed specified nug
    //apply_nugget_build();
    //RChol.copy(R);
    if(ifChooseNug==true)
      nuggetSelectingCholR();
    else
      equationSelectingCholR();
    //equationSelectingCholR();
    if(numRowsR<=numTrend(polyOrder,0)) 
      rcond=0;
    else
      rcond=rcondR;
    //Chol_fact(RChol,chol_info,rcond);
    //if(~(chol_info==0)) assert(chol_info==0); //for debug, do something else for production
    //rcond=rcond_after_Chol_fact(R,RChol);
    //RLU.copy(R);
    //LU_fact(RLU,ipvt_RLU);
    //rcond=rcond_after_LU_fact(R,RLU);
    if(rcond>best_rcond) {
      best_rcond=rcond;
      theta.copy(guess_theta);
      nat_log_corr_len.copy(guess);
    }
  }
  if(constraintType.compare("rcond")==0) { //enforce the condition number bound using rcond
    
    return nat_log_corr_len; //we can't compute analytical derivatives of rcond only the eigenvalues
  }
  else{
    printf("currently the only option for constraint type is \"rcond\"\n");
    assert(constraintType.compare("rcond")==0);
  }
  
  return nat_log_corr_len;
}

/** evaluate the trend function g(x), using specified {Poly, Rot}
Here, g() and x can represent arbitrary (elsewhere represented by
g(x)) or data points (elsewhere represented by G(X)).  The trend
function (unadjusted mean) is dot(g(x),betHat), this returns (matrix)
g(x) for collection of points x. KRD originally implemented this with
a linear trend function */
MtxDbl& KrigingModel::
eval_trend_fn(MtxDbl& g, const MtxInt& poly, 
	      const MtxDbl& rot_or_eul_ang, const MtxDbl& x) const
{
    //static MtxDbl& evalTrendFunc(MtxDbl& g, MtxInt& poly, MtxDbl& rot_or_eul_ang, MtxDbl& x){
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
    return (evaluate_poly_basis(g, poly, xx));
    ///return (LinearRegressionModel::evalBasis(g,poly,rot_or_eul_ang,x));
}


/// evaluate the trend function g(xr), using class members {Poly, Rot}
MtxDbl& KrigingModel::eval_trend_fn(MtxDbl& g, const MtxDbl& xr) const
{
    //MtxDbl& evalTrendFunc(MtxDbl& g, MtxDbl& xr) {
    //printf("KMeTF <");
    eval_trend_fn(g, Poly, Rot, xr);
    //printf("<KMeTF\n");
    return (g);
    
    //return (LinearRegressionModel::evalBasis(g,Poly,Rot,xr));
}

// BMA TODO: These need to be moved to optimizer and then any defauls
// overridden here

void KrigingModel::set_conmin_parameters(OptimizationProblem& opt) const
{
  //set conmin specific parameters for this problem
  //in dot 4.2 the analytical derivative order of objective and constraints must be the same
  if(constraintType.compare("rcond")==0) //rcond constraint can NOT do analytical 1st derivatives
    assert((maxObjDerMode==1)&&(maxConDerMode==1)); //||(maxObjDerMode==3))&&(maxConDerMode==1));
  else {
    printf("currently the only option for constraint type is \"rcond\"\n");
    assert(constraintType.compare("rcond")==0);
  }
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

void KrigingModel::set_direct_parameters(OptimizationProblem& opt) const
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
