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
    numTheta(numVarsr), numPoints(sd.getNPts()), XR(sdBuild.xr), Y(sdBuild.y)
{

  //printf("calling the right KrigingModel constructor\n"); fflush(stdout);

  //if the SurfDataScaler class does what it's supposed to (the only private content in sdBuild that a Model can access are the scaling bits, and only the model can see inside the scaler) the next line will cause an error when you try to compile with it uncommented
  //printf("sd.scaler->mySd.jout=%d\n",scaler.mySd.jout);
  
  // OPTIONS PARSING
  // BMA TODO: max_trials, lower/upper bounds (use toVec helper),
  //           correlation_lengths, optimizationMethod
  ParamMap::const_iterator param_it;

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
  
  // *************************************************************
  // this starts the input section about optimizing or directly
  // scepcifying correlation lengths, it must come after the 
  // scaling section
  // *************************************************************
  
  // current options are none (fixed correl) | sampling (guess) | local | global
  optimizationMethod = "global";
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
  else{ //error checking the input
    cerr << "KrigingModel() unknown optimization_method [" << optimizationMethod << "]  aborting\n";
    assert(0);
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
      if(!(natLogCorrLen(ivarsr)>0.0)) {
	cerr << "For the Kriging Model, correlation lengths must be strictly positive\n.";
	assert(0);
      }

    //printf("unscaled corr lens = [%12.6g",natLogCorrLen(0)); 
    //for(int ivarsr=1; ivarsr<numVarsr; ++ivarsr)
    //printf(", %12.6g",natLogCorrLen(0,ivarsr));
    //printf("]\n");    

    scaler.scaleXrDist(natLogCorrLen); //scale the lengths
    //scaler.scaleXrOther(natLogCorrLen); //error
    //printf("scaled corr lens = [%12.6g",natLogCorrLen(0)); 
    //for(int ivarsr=1; ivarsr<numVarsr; ++ivarsr)
    // printf(", %12.6g",natLogCorrLen(0,ivarsr));
    //printf("]\n");    
    //fflush(stdout);
    
    //compute the natural log of the correlation lengths
    for(int ivarsr=0; ivarsr<numVarsr; ++ivarsr) 
      natLogCorrLen(0,ivarsr)=log(natLogCorrLen(0,ivarsr)); 
    
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
  polyOrder = 1;
  param_it = params.find("order");
  if (param_it != params.end() && param_it->second.size() > 0) {
    polyOrder = std::atoi(param_it->second.c_str()); 
    assert (polyOrder >= 0);
  }
  
  //cout << "order=" << polyOrder << "\n";

  //polyOrder = 2; //for debug
  //main_effects_poly_power(Poly, numVarsr, polyOrder); //for debug
  //commented out for debug
  ifReducedPoly=false;
  param_it = params.find("reduced_polynomial");
  if (param_it != params.end() && param_it->second.size() > 0) 
    if((std::atoi(param_it->second.c_str()))!=0) 
      ifReducedPoly=true;
      
  //cout << "ifReducedPoly=" << ifReducedPoly << "\n";

  if(ifReducedPoly)
    main_effects_poly_power(Poly, numVarsr, polyOrder);
  else
    multi_dim_poly_power(Poly, numVarsr, polyOrder);  

  //check that we have enough data for the selected trend functions
  int needed_points = getNTrend()+1; //+numVarsr;
  if( !(needed_points <= numPoints) ) {
    cerr << "With the selected set of trend functions there are more unknown parameters (" <<  needed_points << ") than there are data points (" << numPoints << ") for the Kriging Model. For the current set of trend functions, you need at least " << needed_points << "data points and twice that is strongly recommended.\n";
    assert(needed_points <= numPoints);
  }

  // *************************************************************
  // this starts the input section HOW to bound the condition 
  // number, this determines which derivatives of the constraint
  // function can be computed analytically so handle that here too
  // *************************************************************
  constraintType="rcond";
  param_it = params.find("constraint_type");
  if (param_it != params.end() && param_it->second.size() > 0)
    if(tolower(param_it->second[0])=='e') //should I use better (whole word) error checking?
      constraintType="eig"; 

  //cout << "contraintType=\"" << constraintType << "\"\n";

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

  assert((0<=num_analytic_obj_ders_in)&&(num_analytic_obj_ders_in<=2)
	 &&(0<=num_analytic_con_ders_in));

  //upper limit to order of analytic derivatives for constraints depends
  //on HOW you want to bound the condition number
  if(constraintType.compare("eig")==0) { //eigenvalue constraint type

    //have implemented analytical first derivatives of eigenvalues
    assert(num_analytic_con_ders_in<=1); 

    numConFunc=numPoints/2; 
    if (numConFunc > 10)
      numConFunc=10;
  }
  else if(constraintType.compare("rcond")==0) { //rcond constraint type    

    //can't do any analytical derivatives (only function values) of rcond
    assert(num_analytic_con_ders_in==0); 
    numConFunc=1;
  }
  else //catch bogus constraint type
    assert((constraintType.compare("eig")==0)||
	   (constraintType.compare("rcond")==0)); 
  
  //convert to the Dakota bitflag convention for derivative orders
  maxObjDerMode=((int) pow(2.0,num_analytic_obj_ders_in+1))-1; //analytical gradients of objective function
  maxConDerMode=((int) pow(2.0,num_analytic_con_ders_in+1))-1; //analytical gradients of constraint function(s)

  //maxCondNum=pow(1024.0,5); 
  maxCondNum=pow(1024.0,4); 
  //maxCondNum=pow(1024.0,5)/32.0;

  // *************************************************************
  // this starts the input section about the nugget which can be
  // used to smooth the data and also decrease the condition 
  // number
  // *************************************************************

  nug=(2*getNTrend()+1.0)/maxCondNum;
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
      cerr << "You can't both auto-select a nugget and use a preset formula" << endl;
      assert(ifChooseNug==false);
    }
    nuggetFormula=std::atoi(param_it->second.c_str()); 
    if(nuggetFormula!=0) {
      switch(nuggetFormula) {
      case 1:
	nug=(2*getNTrend()+1.0)/maxCondNum;
	break;
      case 2:
	nug=2*numPoints/maxCondNum;
	break;
      default:
	cerr << "nugget_formula =" << nuggetFormula << " is not one of the available preset nugget formulas." << endl;
	assert(0);
      }
    }
  }

  param_it = params.find("nugget");
  if (param_it != params.end() && param_it->second.size() > 0) {
    if(!((nuggetFormula==0)&&(ifChooseNug==false))) {
      cerr << "You can do at most 1 of the following (A) auto-select the nugget (minimum needed to satisfy the condition number bound) (B) use one of the preset nugget formulas (C) directly specify a nugget.  The default is not to use a nugget at all (i.e. use a nugget of zero)." << endl;
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
  EulAng.newSize(1, nchoosek(numVarsr, 2)); 
  EulAng.zero();
  gen_rot_mat(Rot, EulAng, numVarsr);
  eval_trend_fn(G, Poly, Rot, XR);
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

  aveDistBetweenPts=pow(numPoints,-1.0/numVarsr);

  //  if(numPoints<=2*numVarsr)
  //  aveDistBetweenPts=1.0;
  //else
  //  aveDistBetweenPts=pow(numPoints-2*numVarsr,-1.0/numVarsr);

  //printf("aveDistBetweenPts=%12.6g\n",aveDistBetweenPts);
  /* Gaussian Process error model has about ~5% confidence (2 std devs away) 
     in what points 4 neighbors away have to say. If points are correlated 
     well at even greater distances then either 
     * that same information will be contained in nearby points OR
     * you shouldn't be using a Gaussian process error model
     KRD */
  double max_corr_length = aveDistBetweenPts*8.0; //*2.0; 

  maxNatLogCorrLen=log(max_corr_length);

  /* Gaussian Process error model has about ~5% confidence (2 std devs) midway
     between neighboring points... i.e. you're 4 std devs away from your 
     nearest neighbor so all sample points are treated as being essentially 
     uncorrelated 
     KRD */
  double min_corr_length = aveDistBetweenPts/4.0; 
  minNatLogCorrLen=log(min_corr_length);
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
    else{
      cerr << "KrigingModel:create() unknown optimization_method [" << optimizationMethod << "]  aborting\n";
      assert(0);
    }
    natLogCorrLen = opt.best_point();
  }

  correlations.newSize(1,numVarsr);
  correlations(0)=0.5*exp(-2.0*natLogCorrLen(0));
  //printf("theta={%g",correlations(0));
  for(int k=1; k<numVarsr; ++k) {
    correlations(k)=0.5*exp(-2.0*natLogCorrLen(k));
    //printf(", %g",correlations(k));
  }

  //printf("}\n");

  //printf("scaled correlations=[%12.6g",correlations(0));
  //for(int ivarsr=1; ivarsr<numVarsr; ++ivarsr)
  //printf(", %12.6g",correlations(0,ivarsr));
  //printf("]\n");

  masterObjectiveAndConstraints(correlations, 1, 0);
  cout << model_summary_string();
  //deallocate matrices we no longer need after emulator has been created

  //temporary variables used by masterObjectiveAndConstraints
  dNbetaHat_dthetaN.clear(); //vector
  temp.clear(); //vector
  temp2.clear(); //vector
  temp3.clear(); //vector
  temp4.clear(); //vector
  temp5.clear(); //vector used to compute hessObj
  d2eps_dthetai_dthetaj.clear(); //vector
  dR_dthetai.clear(); //matrix
  dR_dthetaj.clear(); //matrix used to compute hessObj
  Rinv_dR_dthetai.clear(); //matrix used to compute hessObj
  Rinv_dR_dthetaj.clear(); //matrix used to compute hessObj
  d2R_dthetai_dthetaj.clear(); //matrix used to compute hessObj
  dR_dthetai_Rinv_eps_minus_deps_dthetai.clear(); //matrix used to compute hessObj
  Rinv_dR_dthetai_Rinv_eps_minus_deps_dthetai.clear(); //matrix used to compute hessObj
  Gtran_Rinv_G_inv.clear(); //need this for derivatives of log(det(Gtran_Rinv_G)) but could use it to replace the permanent copy of Gtran_Rinv_G_Chol
  Gtran_Rinv_dR_thetai.clear(); //need this for derivatives of log(det(Gtran_Rinv_G))
  Gtran_Rinv_dR_thetai_Rinv_G.clear(); //need this for derivatives of log(det(Gtran_Rinv_G))

  //variables whose values needed to be retained between sequential call to masterObjectiveAndConstraints for precompute and store strategy to work
  prevObjDerMode=prevConDerMode=0;
  prevTheta.clear(); //row vector 
  allEigVect.clear(); //matrix 
  allEigVal.clear(); //vector
  Z.clear(); //matrix
  R.clear(); //matrix
  G.clear(); //matrix
  deps_dtheta.clear(); //matrix used to compute hessObj
  dR_dtheta_Rinv_eps.clear(); //matrix used to compute hessObj
  destVarianceMLE_dtheta.clear(); //vector
  con.clear(); //vector
  gradObj.clear(); //vector
  gradCon.clear(); //matrix
  hessObj.clear(); //matrix used to compute hessObj
}

std::string KrigingModel::model_summary_string() const {
  MtxDbl temp_out_corr_lengths(1,numVarsr);
  for(int i=0; i<numVarsr; ++i) 
    temp_out_corr_lengths(i)=sqrt(0.5/correlations(i));
  scaler.unScaleXrDist(temp_out_corr_lengths);
  
  std::ostringstream oss;
  oss << "KM: #pts=" << numPoints <<"; Correlation lengths=(" << temp_out_corr_lengths(0);
  for(int i=1; i<numVarsr; ++i)
    oss << ", " << temp_out_corr_lengths(i);
  oss << "); unadjusted variance=" << estVarianceMLE * scaler.unScaleFactorVarY() << "; log(likelihood)=" << likelihood << "; the trend is a ";
  if(polyOrder>1) {
    if(ifReducedPoly==true)
      oss << "reduced_";
    else oss <<"full ";
  }
  oss << "polynomial of order=" << polyOrder << 
    "; rcond(R)=" << rcondR << "; rcond(Gtran_Rinv_G)=" << rcondGtran_Rinv_G 
      << "; nugget=" << nug << ".\n";
	
  oss << "Beta= (" << betaHat(0);
  int ntrend=getNTrend(); 
  for(int i=1; i<ntrend; ++i)
    oss << "," << betaHat(i);
  oss << ")\n";

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
  int ntrend = getNTrend(); 
  MtxDbl g(1, ntrend), r(1, numPoints);

  /*
  printf("double evaluate()\n");
  printf("xr=[%20.14g", xr(0));
  for(int i=1; i<numVarsr; ++i)
    printf(", %20.14g",xr(i));
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
      
  if(nug>0.0)
    apply_nugget_eval(r);
  
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
      y(i)=singular_y;
    return y;
  }
  
  //assert( (numVarsr == xr.getNCols()) && (xr.getNRows() == 1) );
  int ntrend = getNTrend(); 
  MtxDbl g(nrowsxr, ntrend), r(nrowsxr, numPoints);

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
    printf("] Y(%3d)=%12.6g\n",i,Y(i));
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

  int ntrend=getNTrend();

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
  
  MtxDbl r(nrowsxr,numPoints);
  correlation_matrix(r, xr_scaled);
  apply_nugget_eval(r);
  MtxDbl d1r(nrowsxr,numPoints);
  MtxDbl temp_vec(nrowsxr);


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
      d1y(ipt,ider)=(d1y(ipt,ider)+temp_vec(ipt))*d1y_unscale_factor;
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

  MtxDbl r(nrowsxr,numPoints);
  correlation_matrix(r, xr);
  apply_nugget_eval(r);
  MtxDbl d1r(nrowsxr,numPoints);
  MtxDbl d2r(nrowsxr,numPoints);
  MtxDbl temp_vec(nrowsxr);

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
      d2y(ipt,ider)=(d2y(ipt,ider)+temp_vec(ipt))*d2y_unscale_factor;
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
  int ntrend = getNTrend(); 
  MtxDbl g_minus_r_Rinv_G(1, ntrend), r(1, numPoints);

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

  MtxDbl tempa(numPoints);
  MtxDbl tempb(ntrend);

  if(nug>0.0) 
    apply_nugget_eval(r);

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
  int ntrend = getNTrend(); 
  MtxDbl g_minus_r_Rinv_G(nrowsxr, ntrend), r(nrowsxr, numPoints);

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
  MtxDbl tempa(numPoints,nrowsxr); 
  MtxDbl tempb(ntrend,nrowsxr);
  double var_unscale_factor=scaler.unScaleFactorVarY();

  if(nug>0.0)
    apply_nugget_eval(r);
  
  solve_after_Chol_fact(tempa,RChol,r,'T');

  matrix_mult(g_minus_r_Rinv_G,r,Rinv_G,1.0,-1.0);
  solve_after_Chol_fact(tempb,Gtran_Rinv_G_Chol,g_minus_r_Rinv_G,'T');

  for(i=0; i<nrowsxr; ++i) {
    //saved 2*nrowsxr loops
    adj_var(i)=1.0-r(i,0)*tempa(0,i)+g_minus_r_Rinv_G(i,0)*tempb(0,i);

    for(j=1; j<numPoints; ++j)
      adj_var(i)-=r(i,j)*tempa(j,i); //looks a lot like matrix mult but only N^2 ops

    for(j=1; j<ntrend; ++j)
      adj_var(i)+=g_minus_r_Rinv_G(i,j)*tempb(j,i); //looks a lot like matrix mult but only N^2 ops

    adj_var(i)*=estVarianceMLE*var_unscale_factor;
  }

  for(i=0; i<nrowsxr; ++i)
    if(adj_var(i)<0.0)
      adj_var(i)=0.0;
    else if(!(adj_var(i)>=0.0))
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


/** this function scales the input matrix, r, in place by 1.0/(1.0+nug), 
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

  int nelem=r.getNElems();
  double temp_dbl=1.0/(1.0+nug);

  for(int ij=0; ij<nelem; ++ij) 
    r(ij)*=temp_dbl;

  //printf("apply_nugget_eval temp=%g\n",temp);

  return r;
}

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
  
  int nrowsR=R.getNRows();
  //assert(nrowsR==R.getNCols());
  int nelemsR=nrowsR*nrowsR;

  double temp_dbl=1.0/(1.0+nug);
  int ij;
  for(ij=0; ij<nelemsR; ++ij)
    R(ij)*=temp_dbl;
  
  //the "paranoid" part of my mind wonders if there would be less round off
  //error if I added the nugget to the diagonal BEFORE scaling, the pragmatic
  //part of my mind says it shouldn't matter and doing it this way is faster
  temp_dbl*=nug;
  for(ij=0; ij<nrowsR; ++ij)
    R(ij,ij)+=temp_dbl; 
    
  return;
}

/** r (lower case r) is the correlation matrix between the
    interpolation points and data points, it used to EVALUATE but not
    construct the emulator's Gaussian process error model
    i.e. E(y(xr)|Y(XR))=g(xr)*betaHat+r*R^-1*eps where
    eps=(Y-G(XR)*betaHat), to be more specific
    r(i,j)=r(xr(i),XR(j))=exp(sum_k -theta(k)*(xr(i,k)-XR(j,k))^2) KRD
    wrote this */
MtxDbl& KrigingModel::correlation_matrix(MtxDbl& r, const MtxDbl& xr) const
{
  int nrowsXR=XR.getNRows(); //data points used to build model
  int nrowsxr=xr.getNRows(); //points at which we are evalutating the model
  //assert(xr.getNCols()==numVarsr);
  r.newSize(nrowsxr,nrowsXR);
  int nelemr=nrowsxr*nrowsXR;
  int i,j,k, ij;
  double temp_double;

  //REVACC_F77(r.ptr(0),xr.ptr(0),XR.ptr(0),correlations.ptr(0),
  //     &nrowsxr,&nrowsXR,&numVarsr);

  /*
  printf("**********************************************************\n");
  printf("**********************************************************\n");
  for(i=0; i<nrowsXR; i++) {
    printf("XR(%d,:)={",i);
    for(j=0; j<numVarsr; j++) printf(" %g",XR(i,j));
    printf(" }; Y(%d)=%g\n",i,Y(i));
  }
  printf("**********************************************************\n");
  printf("**********************************************************\n");
  */
  
  if(numVarsr==1) {
    //optimized for when there is only 1 output variable
    for(j=0; j<nrowsXR; j++)
      for(i=0; i<nrowsxr; i++){
	temp_double=xr(i,0)-XR(j,0);
	r(i,j)=exp(-correlations(0)*temp_double*temp_double); 
    }
  } else if(nrowsxr==1) {
    //"optimized" for when there is only 1 evaluation point (Save loops, save dereferences)
    
    //k=0 was pulled out from below to avoid doing an extra loop just to 
    //initialize all of r to zero
    for(j=0; j<nrowsXR; j++) {
      temp_double=xr(0)-XR(j,0);
      r(j)=-correlations(0)*temp_double*temp_double; //=- is correct
    }
  
    //all but first and last k
    for(k=1;k<numVarsr-1;k++)
      for(j=0; j<nrowsXR; j++) {
	temp_double=xr(k)-XR(j,k);
	r(j)-=correlations(k)*temp_double*temp_double; //-= is correct
      }
  
    //this value of k was pulled out of above to save doing an extra loop 
    //for just the exp() operation
    k=numVarsr-1; 
    for(j=0; j<nrowsXR; j++) {
      temp_double=xr(k)-XR(j,k);
      r(j)=exp(r(j)-correlations(k)*temp_double*temp_double); 
    }
  } else {
    //"optimized" for when there is more than 1 output variable
    
    //k=0 was pulled out from below to avoid doing an extra loop just to 
    //initialize all of r to zero
    for(j=0; j<nrowsXR; j++)
      for(i=0; i<nrowsxr; i++){
	temp_double=xr(i,0)-XR(j,0);
	r(i,j)=-correlations(0)*temp_double*temp_double; //=- is correct
      }
  
    //all but first and last k
    for(k=1;k<numVarsr-1;k++)
      for(j=0; j<nrowsXR; j++)
	for(i=0; i<nrowsxr; i++){
	  temp_double=xr(i,k)-XR(j,k);
	  r(i,j)-=correlations(k)*temp_double*temp_double; //-= is correct
	}
  
    //this value of k was pulled out of above to save doing an extra loop 
    //for just the exp() operation
    k=numVarsr-1; 
    for(j=0; j<nrowsXR; j++)
      for(i=0; i<nrowsxr; i++){
	temp_double=xr(i,k)-XR(j,k);
	r(i,j)=exp(r(i,j)-correlations(k)*temp_double*temp_double); 
      }
  }
  return r;
}

///k is the variable/dimension not the point
MtxDbl& KrigingModel::dcorrelation_matrix_dxI(MtxDbl& dr, const MtxDbl& r, 
					      const MtxDbl& xr, int Ider) const
{
  int nrowsxr=xr.getNRows();
  assert((r.getNRows()==nrowsxr)&&(r.getNCols()==numPoints)&&
	 (xr.getNCols()==numVarsr)&&(0<=Ider)&&(Ider<numVarsr));
  dr.newSize(nrowsxr,numPoints);

  double temp_dbl=-2.0*correlations(Ider);
  for(int jpt=0; jpt<numPoints; ++jpt)
    for(int ipt=0; ipt<nrowsxr; ++ipt)
      dr(ipt,jpt)=temp_dbl*r(ipt,jpt)*(xr(ipt,Ider)-XR(jpt,Ider));

  return dr;
}

MtxDbl& KrigingModel::d2correlation_matrix_dxIdxK(MtxDbl& d2r, const MtxDbl& drI, const MtxDbl& r, const MtxDbl& xr, int Ider, int Kder) const
{
  int nrowsXR=XR.getNRows(); //data points used to build model
  int nrowsxr=xr.getNRows(); //points at which we are evalutating the model
  d2r.newSize(nrowsxr,nrowsXR);
///k is the variable/dimension not the point

  assert((r.getNRows()==nrowsxr)&&(r.getNCols()==numPoints)&&
	 (xr.getNCols()==numVarsr)&&(0<=Kder)&&(Kder<numVarsr));

  double neg_two_theta_K=-2.0*correlations(Kder);
  if(Ider==Kder) 
    for(int jpt=0; jpt<numPoints; ++jpt)
      for(int ipt=0; ipt<nrowsxr; ++ipt)
	d2r(ipt,jpt)=neg_two_theta_K*((xr(ipt,Kder)-XR(jpt,Kder))*drI(ipt,jpt)+r(ipt,jpt));
  else
    for(int jpt=0; jpt<numPoints; ++jpt)
      for(int ipt=0; ipt<nrowsxr; ++ipt)
	d2r(ipt,jpt)=neg_two_theta_K*(xr(ipt,Kder)-XR(jpt,Kder))*drI(ipt,jpt);

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
  R.newSize(numPoints,numPoints);
  double *R_ptr=R.ptr(0);

  double Rij_temp;
  int ij=0;
  for(int j=0; j<numPoints-1; ++j) {
    R(j,j)=1.0;
    for(int i=j+1; i<numPoints; ++i) {
      Rij_temp=exp(Ztheta(ij));
      R(i,j)=Rij_temp;
      R(j,i)=Rij_temp;
      ++ij;
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
  register int ijk=0;
  double *Z_ptr=Z.ptr(0); //done for speed
  const double *XR_k_ptr; //done for speed
  for(int k=0; k<ncolsXR; k++) {
    XR_k_ptr=XR.ptr(0,k);
    for(int j=0; j<nrowsXR-1; j++)
      for(int i=j+1; i<nrowsXR; i++) {
	mult_term=XR_k_ptr[i]-XR_k_ptr[j];
	Z_ptr[ijk++]=-mult_term*mult_term;
      }
  }
  
  return Z;
}

/** Warning it's up to the person calling this function to ensure that
    R is the correlation Matrix generated as R=exp(Z*theta) where
    theta is the correlation vector and Z(ij,k)=-(XR(i,k)-XR(j,k))^2. I
    (KRD) did it this way because it is faster and the only place it
    is currently being used is in void
    KrigingProblem::objectiveAndGradient().  the Z matrix is produced
    by KrigingModel::gen_Z_matrix() and R is produced by MtxDbl
    KrigingModel::corrMtx(const int nrowsXR, const MtxDbl& Z, const
    VecDbl& theta) KRD wrote this. Z is a member variables so it doesn't 
    need to be passed in, R needs to be passed in so we can evaluate 
    second derivatives with the same function... dR_dthetak=-Z(:,k).*R(:) 
    (MATLAB notation), and d2R_dthetai_dthetaj=-Z(:,j).*dR_dthetai;*/
MtxDbl& KrigingModel::dcorrMtx_dthetak(MtxDbl& dR_dthetak, const MtxDbl& R_local, 
				       const int k)
{
  dR_dthetak.newSize(numPoints,numPoints);
  double dRij_dthetak_temp;
  int nrowsZ=Z.getNRows();
  const double *Zk_ptr=Z.ptr(0,k); //do this to avoid one extra dereference per element of R_local
  int ij=0;
  for(int j=0; j<numPoints-1; ++j) {
    dR_dthetak(j,j)=0.0;
    for(int i=j+1; i<numPoints; ++i) {
      dRij_dthetak_temp=Zk_ptr[ij]*R_local(i,j);
      dR_dthetak(i,j)=dRij_dthetak_temp;
      dR_dthetak(j,i)=dRij_dthetak_temp;
      ++ij;
    }      
  }
  dR_dthetak(numPoints-1,numPoints-1)=0.0;

  return dR_dthetak;
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
  //member variable: MtxDbl dNbetaHat_dthetaN(ntrend)
  //member variable: MtxDbl temp(ntrend)
  //member variable: MtxDbl temp2(numPoints)
  //member variable: MtxDbl temp3(numPoints)
  //member variable: MtxDbl temp4(numPoints)
  //member variable: MtxDbl temp5(numPoints)
  //member variable: MtxDbl d2eps_dthetai_dthetaj(numPoints)
  //member variable: MtxDbl dR_dthetai(numPoints,numPoints)
  //member variable: MtxDbl dR_dthetaj(numPoints,numPoints)
  //member variable: MtxDbl Rinv_dR_dthetai(numPoints,numPoints)
  //member variable: MtxDbl Rinv_dR_dthetaj(numPoints,numPoints)  
  //member variable: MtxDbl d2R_dthetai_dthetaj(numPoints,numPoints)
  //member variable: MtxDbl dR_dthetai_Rinv_eps_minus_deps_dthetai(numPoints,numTheta)
  //member variable: MtxDbl Rinv_dR_dthetai_Rinv_eps_minus_deps_dthetai(numPoints,numTheta)
  


  //these are private and their values need to be retained between sequential calls to masterObjectiveAndConstraints, that is no other function (other than the create) can access them, these get deallocated at the end of create
  //member variable: int prevObjDerMode
  //member variable: int prevConDerMode
  //member variable: MtxDbl prevTheta(1,numTheta)
  //member variable: MtxDbl allEigVect(numPoints,numPoints)
  //member variable: MtxDbl allEigVal(numPoints)
  //member variable: MtxDbl Z(numPoints*numPoints,numTheta)
  //member variable: MtxDbl R(numPoints,numPoints)
  //member variable: MtxDbl G(numPoints,ntrend)
  //member variable: MtxDbl deps_dtheta(numPoints,numTheta)
  //member variable: MtxDbl dR_dtheta_Rinv_eps(numPoints,numTheta)
  //member variable: MtxDbl destVarianceMLE_dtheta(numTheta)


  //keep these around after emulator creation so we can evaluate the emulator
  //member variable: MtxDbl rhs(numPoints)
  //member variable: MtxDbl RLU(numPoints,numPoints)
  //member variable: MtxInt ipvt_RLU
  //member variable: MtxDbl Rinv(numPoints,numPoints) //needed to eval integral of adjusted variance, capability not yet added
  //member variable: MtxDbl Rinv_G(numPoints,ntrend)
  //member variable: MtxDbl Gtran_Rinv_G_LU(ntrend,ntrend)
  //member variable: MtxInt ipvt_Gtran_Rinv_G_LU(ntrend)

  //if theta was the same as the last time we called this function than we can reuse some of the things we calculated last time
  
  int i, j;

  if(prevTheta.getNElems()!=numTheta) {
    prevTheta.newSize(1,numTheta);
    prevObjDerMode=prevConDerMode=0; 
  }
  else
    for(i=0; i<numTheta; ++i) 
      if(prevTheta(i)!=theta(i)) {
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
      prevTheta(i)=theta(i); 


  int k;
  int ntrend=Poly.getNRows();
  int chol_info;

  if((prevObjDerMode==0)&&(prevConDerMode==0)) {
    R.newSize(numPoints,numPoints);
    correlation_matrix(theta); //fills member variable R as exp(Z*theta) where Z is a member variable
    apply_nugget_build(); //modify R by nug in place
  }

  if((prevObjDerMode==0)&&((1<=obj_der_mode)||(1<=con_der_mode))) {


    //perform LU decomposition of R and calculate the determinant of R, replaced LU with Cholesky
    //http://en.wikipedia.org/wiki/Determinant#Determinant_from_LU_decomposition
    RChol.copy(R);
    chol_info=0;
    Chol_fact(RChol,chol_info,rcondR); //preconditioned Cholesky, when Kriging is not gradient enhaced R won't need preconditiong since it has all 1's on the diagonals
    //printf("Chol\n");
    //assert(chol_info==0);  //here for debug, decide what to do about it later

    double log_determinant_R = 0.0; //need to do this to avoid underflow error for large numbers of points, log(0)=-inf
    for (int ipt = 0; ipt < numPoints; ++ipt) 
      log_determinant_R += log(RChol(ipt,ipt)); 
    log_determinant_R *= 2.0; //only for Cholesky factorization of R 

    //determinant_R=fabs(determinant_R); //KRD added fabs because "The determinant of a positive definite matrix is always positive" http://mathworld.wolfram.com/PositiveDefiniteMatrix.html and det(R)=det(pivot Mtx)*det(L)*det(U); det(L)=1, det(U) is what we calculated above and det(pivot Mtx)=+/- 1

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
    Rinv_G.newSize(numPoints,ntrend); //precompute and store
    solve_after_Chol_fact(Rinv_G,RChol,G);
    //solve_after_LU_fact(Rinv_G,RLU,ipvt_RLU,G,'N','N'); //O(N^3) ops
    Gtran_Rinv_G_Chol.newSize(ntrend,ntrend);
    matrix_mult(Gtran_Rinv_G_Chol,G,Rinv_G,0.0,1.0,'T','N');
    //double rcondGtran_Rinv_G;
    Chol_fact(Gtran_Rinv_G_Chol,chol_info,rcondGtran_Rinv_G);

    double log_determinant_Gtran_Rinv_G=0.0;
    for (int itrend = 0; itrend < ntrend; ++itrend)
      log_determinant_Gtran_Rinv_G += log(Gtran_Rinv_G_Chol(itrend,itrend)); 
    log_determinant_Gtran_Rinv_G *= 2.0; //only for Cholesky factorization of R 
    //if(~(chol_info==0)) assert(chol_info==0);  //for debug, do something else for production

    temp.newSize(ntrend);
    matrix_mult(temp, Rinv_G, Y, 0.0, 1.0, 'T', 'N');
    betaHat.newSize(ntrend);
    solve_after_Chol_fact(betaHat,Gtran_Rinv_G_Chol,temp); //O(ntrend^2) ops

    temp2.copy(Y); //this will be eps=epsilon=Y-G(XR)*betaHat, but use 
    //variable temp2 because we would only need the variable "eps" for these 
    //5 lines of code (not counting comments) and we want to save space, 
    //afterwards we will only need R^-1*eps which is stored in "rhs"
    matrix_mult(temp2, G, betaHat, 1.0, -1.0, 'N', 'N'); //eps=Y-G(XR)*betaHat
    rhs.newSize(numPoints);
    solve_after_Chol_fact(rhs,RChol,temp2);




    //it's actually the log likelihood, which we want to maximize
    //likelihood = -0.5*(numPoints*(log(4.0*acos(0.0))+log(estVarianceMLE)+1)
    //		       +log(determinant_R)); //from Koehler and Owen 

#ifdef __NKM_UNBIASED_LIKE__
    //derived following: C. E. Rasmussen & C. K. I. Williams, Gaussian Processes for Machine Learning, the MIT Press, 2006, ISBN 026218253X. c 2006 Massachusetts Institute of Technology. www.GaussianProcess.org/gpml...  we assume a "vague prior" (i.e. that we don't know anything) for betaHat, then like "Koehler and Owen" we replace the covariance matrix K with (unadjusted variance)*R (where R is the correlation matrix) and find unadjusted variance and betaHat through maximum likelihood.

    //the unbiased estimate of unadjusted variance
    estVarianceMLE = dot_product(temp2,rhs)/(numPoints-ntrend); 

    //the "per point" unbiased log(likelihood)
    likelihood = -0.5*(log(estVarianceMLE)+(log_determinant_R+log_determinant_Gtran_Rinv_G)/(numPoints-ntrend)); 
#else
    //derived the "Koehler and Owen" way (assumes we know the trend function, and is therefore biased, but usally seems to work better for surrogate based optimization)

    //the estimate of unadjusted variance
    estVarianceMLE = dot_product(temp2,rhs)/numPoints; //the "Koehler and Owen" way

    //the "per point" log(likelihood)
    likelihood = -0.5*(log(estVarianceMLE)+log_determinant_R/numPoints); 
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
    con.newSize(numConFunc);
    
    if(constraintType.compare("rcond")==0) { //use rcond (and maybe its numerical derivatives) to bound the 
      //condition number
      
      assert((1<=prevObjDerMode)&&(numConFunc==1)); //make sure we have calculated rcondR already
      //con(0)=1.0-rcondR*3.0*maxCondNum;  //have seen rcond as low as about 1/3 of the true value
      con(0)=1.0-rcondR*maxCondNum;  //have seen rcond as low as about 1/3 of the true value
    }
    else if(constraintType.compare("eig")==0) { //use eigenvalues (and maybe their analytical derivatives) to
      //bound the condition number
      allEigVect.newSize(numPoints,numPoints);
      allEigVal.newSize(numPoints);

      eig_sym(allEigVect, allEigVal, R); //R is correct
      for(int icon=0; icon<numConFunc; icon++)
	con(icon)=(allEigVal(numPoints-1)+nug)/maxCondNum-(allEigVal(icon)+nug);
        //con(icon)=1.0-maxCondNum*(allEigVal(icon)+nug)/(allEigVal(numPoints-1)+nug);
        //con(icon)=1.0/maxCondNum-(allEigVal(icon)+nug);
    }
    else
      assert((constraintType.compare("eig")==0)||(constraintType.compare("rcond")==0));
    
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
    //if you want gradients of the objective AND constraint functions 
    //do them together below, so they can share the calculation of 
    //dR_thethai.  However if you want gradients of the contraint 
    //functions BUT DO NOT WANT the gradient of the objective function 
    //THEN calculate the gradients of the constraints here
    assert(constraintType.compare("eig")==0); //we can compute analytical 
    //derivatives of only the eigenvalues (we can't compute analytical
    //derivatives of rcondR)

    gradCon.newSize(numConFunc,numTheta);
    dR_dthetai.newSize(numPoints,numPoints);
    temp2.newSize(numPoints);
    temp3.newSize(numPoints);
    temp4.newSize(numPoints); //FOR CALCULATION OF THE GRADIENTS OF THE 
    //CONSTRAINT FUNCTIONS this holds the eigenvector corresponding to the
    //the largest eigenvalue    
    allEigVect.getCols(temp4,numPoints-1); 

    double temp_con;
    //caclulate the gradient of the constraints...
    //constraint(i)=(largest_eigenvalue/Condition_number_bound - 
    //               the_ith_smallest_eigenvalue) < 0 
    //obviously if this holds true for the first smallest eigenvalue 
    //it will hold true for the 2nd, 3rd, etc smallest eigenvalues 
    //but I use i=1->numConFunc in case any of the small eigenvalues 
    //trade positions (note that C arrays start at 0 instead of 1)
    //where numConFunc=min(numPoints/2,10) (a heuristic).  Of course we 
    //need the analytical derivatives of the eigenvalues.  The formula 
    //for exact arbitray order derivatives is found in 
    //  Jankovic, Exact nth Derivatives of Eigenvalues and Eigenvectors.
    //     Journal of Guidance, Control, and Dynamics, Vol 17. No 1,
    //     January-February 1994 pp 136-144
    //I (KRD) found this paper to be difficult to read but highly useful 
    //once it is deciphered.  To help with the deciphering, Eqns {1, 
    //20, and 22->24} are the only ones you will typically need for 
    //_linear_ algebra (example 1 is useful too), if your eigenvectors 
    //are normalized, then the denominator in equation 24 is -1 for 
    //linear algebra
    double dlarge_eigval_dthetai, dsmall_eigval_dthetai;
    for(i=0; i<numTheta; ++i) {
      dcorrMtx_dthetak(dR_dthetai, R, i);
      apply_nugget_eval(dR_dthetai); //apply nugget EVAL is correct, the additive term is constant with respect to theta so drops out

      dlarge_eigval_dthetai=
	dot_product(temp4,matrix_mult(temp2,dR_dthetai,temp4,0.0,1.0));
      
      temp_con=	dot_product(temp4,
      	    matrix_mult(temp2,dR_dthetai,temp4,0.0,1.0))/
      maxCondNum;
      for(int icon=0; icon<numConFunc; icon++) {
	allEigVect.getCols(temp3,icon);
	//dsmall_eigval_dthetai=
	//dot_product(temp3,matrix_mult(temp2,dR_dthetai,temp3,0.0,1.0));
	//gradCon(icon,i)=-dsmall_eigval_dthetai;
	//gradCon(icon,i)=-maxCondNum*
	//(dsmall_eigval_dthetai*(allEigVal(numPoints-1)+nug)-
	// dlarge_eigval_dthetai*(allEigVal(icon)+nug))/
	//((allEigVal(numPoints-1)+nug)*(allEigVal(numPoints-1)+nug));
		gradCon(icon,i)=temp_con-
	dot_product(temp3,
	      matrix_mult(temp2,dR_dthetai,temp3,0.0,1.0));
      }
    }

    prevConDerMode=2; //increase prevConDerMode to current value
    if((obj_der_mode<=prevObjDerMode)&&(con_der_mode==2)) {
      //early exit because we have everything we need
      return;
    }
  }

  //******************************************************************
  //******************************************************************
  // calculating the gradient of the objective function (and optionally
  // the constraint functions) starts here
  //******************************************************************
  //******************************************************************
  if((prevObjDerMode<2)&&(2<=obj_der_mode)) {
    gradObj.newSize(numTheta); 

    if(2<=con_der_mode) 
      gradCon.newSize(numConFunc,numTheta);

    dR_dthetai.newSize(numPoints,numPoints);

    destVarianceMLE_dtheta.newSize(numTheta);
    //since we only need the TRACE of the product of R^-1*dR_dthetai we 
    //can reshape R^-1 and dR_dtheta into vectors and take their dot 
    //product consting a N^3+numTheta*N^2 ops (this N^3 comes from 
    //inverting R) instead of numTheta*N^3 ops (these N^3 comes from 
    //LU_solve(R\dR_dthetai) for each i=0->numTheta)    
    //Rinv.copy(RLU); //Rinv=R^-1, NOT a pseudo inverse
    //inverse_after_LU_fact(Rinv, ipvt_RLU); //the only O(N^3) operarion
    Rinv.copy(RChol);
    inverse_after_Chol_fact(Rinv);
    Gtran_Rinv_G_inv.copy(Gtran_Rinv_G_Chol);
    inverse_after_Chol_fact(Gtran_Rinv_G_inv);
    //needed to calculate the gradients of the objective function
    
    assert(ntrend>0); //for debug
    temp.newSize(ntrend); //FOR CALCULATION OF THE GRADIENT OF THE OBJECTIVE 
    //FUNCTION this holds -G^T*R^-1*dR_dthetai*R^-1*(Y-G*betaHat)
    //for calculation of the hessian of the objective function it holds 
    //something different.  To avoid the cost of constantly allocating and 
    //deallocating memory, it should ALWAYS be a vector with numPoints elements

    temp2.newSize(numPoints); //FOR CALCULATION OF THE GRADIENT OF THE 
    //CONSTRAINT functions this holds dR_dthetai times the largest eigenvector
    //FOR CALCULATION OF THE GRADIENT OF THE OBJECTIVE FUNCTION this holds 
    //dR_dthetai*R^-1*eps where eps=Y-G*betaHat. For calculation of the 
    //hessian of the objective function this holds something different.  To 
    //avoid the cost of constantly allocating and deallocating memory, it 
    //should ALWAYS be a vector with numPoints elements
    
    temp3.newSize(numPoints); //FOR CALCULATION OF THE GRADIENT OF THE 
    //CONSTRAINT FUNCTION this holds the eigen vector corresponding to
    //a small eigenvalue (ratio of largest to smallest eigenvalue is 
    //constrained but small eigenvalues can swap places) FOR CALCULATION OF 
    //THE GRADIENT OF THE OBJECTIVE FUNCTION temp3 holds a column of 
    //deps_dtheta. FOR CALCULATION OF THE HESSIAN OF THE OBJECTIVE FUNCTION 
    //temp3 holds something else.  To avoid the cost of constantly allocating 
    //and deallocating memory, it should ALWAYS be a vector with numPoints 
    //elements

    if((prevConDerMode<2)&&(2<=con_der_mode)) {
      temp4.newSize(numPoints); //FOR CALCULATION OF THE GRADIENTS OF THE 
      //CONSTRAINT FUNCTIONS this holds the eigenvector corresponding to the
      //the largest eigenvalue
    }

    dNbetaHat_dthetaN.newSize(ntrend); //This is another temporary variable 
    //(we do not need to retain its state between sequential calls to this 
    //function) but I gave it a name reflecting its contents. FOR CALCULATION 
    //OF THE GRADIENT OF THE OBJECTIVE FUNCTION dNbetaHat_dthetaN holds the 
    //first derivative of the vector of least squares coefficients with 
    //respect to theta(i). FOR CALCULATION OF THE HESSIAN OF THE OBJECTIVE 
    //FUNCTION it holds the mixed second derivative of the least squares 
    //coefficients with respect to theta(i) and theta(j).  

    dR_dthetai.newSize(numPoints,numPoints); //you can consider this a temp 
    //matrix, newSize()-ing is either cheap (if it's already the right size)
    //or necessary, we declared this MtxDbl to be a member variable because
    //we will often needed it and that let's us avoid constantly allocating 
    //and deallocating memory.

    if(4<=maxObjDerMode) {
      //if we are eventually going to calculate hessian of the objective 
      //function then I want to save some gradient information for use in
      //calculating the hessian, and I need to allocate space for it
      dR_dtheta_Rinv_eps.newSize(numPoints,numTheta);
      deps_dtheta.newSize(numPoints,numTheta);
    }

    double trace_Rinv_dR_dthetai;
    double trace_Gtran_Rinv_G_inv_Gtran_Rinv_dR_dthetai_Rinv_G;
    double temp_con;
    double dlarge_eigval_dthetai, dsmall_eigval_dthetai;
    if((prevConDerMode<2)&&(2<=con_der_mode)) {
      //only do this if we also need the gradients of the constraint functions
      allEigVect.getCols(temp4,numPoints-1); //temp4 is the eigenvector with 
      //the largest eigenvalue, hereafter refered to as the largest eigenvector
    }

    //use the precomputed and stored quantities to quickly calculate the 
    //gradient of the objective and constraint functions
    for (i = 0; i<numTheta; i++) {
      dcorrMtx_dthetak(dR_dthetai, R, i); //dR_dthetai and R are correct, at this point dR_dthetai is dR_dthetai
      apply_nugget_eval(dR_dthetai); //apply nugget EVAL is correct, the additive term is constant with respect to theta so drops out
      
      if((prevConDerMode<2)&&(2<=con_der_mode)) {
	assert(constraintType.compare("eig")==0); //we can compute analytical 
	//derivatives of only the eigenvalues (we can't compute analytical
	//derivatives of rcondR)

	//caclulate the gradient of the constraints...
	//constraint(i)=(largest_eigenvalue/Condition_number_bound - 
	//               the_ith_smallest_eigenvalue) < 0 
	//obviously if this holds true for the first smallest eigenvalue 
	//it will hold true for the 2nd, 3rd, etc smallest eigenvalues 
	//but I use i=1->numConFunc in case any of the small eigenvalues 
	//trade positions (note that C arrays start at 0 instead of 1)
	//where numConFunc=min(numPoints/2,10) (a heuristic).  Of course we 
	//need the analytical derivatives of the eigenvalues.  The formula 
	//for exact arbitray order derivatives is found in 
	//  Jankovic, Exact nth Derivatives of Eigenvalues and Eigenvectors.
	//     Journal of Guidance, Control, and Dynamics, Vol 17. No 1,
	//     January-February 1994 pp 136-144
	//I (KRD) found this paper to be difficult to read but highly useful 
	//once it is deciphered.  To help with the deciphering, Eqns {1, 
	//20, and 22->24} are the only ones you will typically need for 
	//_linear_ algebra (example 1 is useful too), if your eigenvectors 
	//are normalized, then the denominator in equation 24 is -1 for 
	//linear algebra

	//dlarge_eigval_dthetai=
	//dot_product(temp4,matrix_mult(temp2,dR_dthetai,temp4,0.0,1.0));
	//	for(int icon=0; icon<numConFunc; icon++) {
	//allEigVect.getCols(temp3,icon);
	//dsmall_eigval_dthetai=
	//  dot_product(temp3,matrix_mult(temp2,dR_dthetai,temp3,0.0,1.0));
	//gradCon(icon,i)=-dsmall_eigval_dthetai;
	  //gradCon(icon,i)=-maxCondNum*
	  //(dsmall_eigval_dthetai*(allEigVal(numPoints-1)+nug)-
	  // dlarge_eigval_dthetai*(allEigVal(icon)+nug))/
	  //((allEigVal(numPoints-1)+nug)*(allEigVal(numPoints-1)+nug));
	//}

	temp_con=
	dot_product(temp4,
	      matrix_mult(temp2,dR_dthetai,temp4,0.0,1.0))/
	maxCondNum;
	for(int icon=0; icon<numConFunc; icon++) {
	allEigVect.getCols(temp3,icon); //temp3 is a "small eigen vector"
	gradCon(icon,i)=temp_con-
	  dot_product(temp3,
		matrix_mult(temp2,dR_dthetai,temp3,0.0,1.0));
	}
      }

      //  the math
      //  R=(1+nug)^-1*(R+nug*I), typically nug=0
      //  betaHat=(G^T*R^-1*G)^-1*(G^T*R^-1*Y)
      //  d(A^-1)_dthetai=-A^-1*dA_dthetai*A^-1 (where A is either R or G^T*R^-1*G)
      //  dbetaHat_dthetai=-(G^T*R^-1*G)^-1*G^T*R^-1*dR_dthetai*R^-1*eps
      //  eps=Y-G*betaHat
      //  deps_dthetai=-G*dbetaHat_dthetai
      //  estVarMLE=1/N*eps^T*R^-1*eps
      //  destVarMLE_dthetai=1/N*[2*deps_dthetai^T*R^-1*eps-eps^T*R^-1*dR_dthetai*R^-1*eps]
      //  ddetR_dthetai=det(R)*Trace[R^-1*dR_dthetai]
      //  dlogdetR_dthetai=Trace[R^-1*dR_dthetai] (store R^-1 so can evaluate this in O(N^2) instead of O(N^3) ops)
      //  obj=0.5*[N*(log(2*pi)+log(estVarMLE)+1)+log(det(R))]
      //  dobj_dthetai=0.5*[N*destVarMLE_dthetai/estVarMLE+dlogdetR_dthetai]
      //
      //  we have already precomputed and stored
      //  G^T*R^-1
      //  (G^T*R^-1*G)^-1
      //  R^-1*G*betaHat-R^-1*Y
      //  R^-1*eps
      //  R^-1
      //
      //  will loop over i to fill up the gradient of the objective function

      //dR_dthetai*R^-1*(Y-G*betaHat)
      matrix_mult(temp2, dR_dthetai, rhs, 0.0, 1.0); 
      if(4<=maxObjDerMode) {
	//only store this if we're going to calculate the objective 
	//function's hessian
	dR_dtheta_Rinv_eps.putCols(temp2,i);
      }
      
      //-G^T*R^-1*dR_dthetai*R^-1*(Y-G*betaHat)
      matrix_mult(temp, Rinv_G, temp2, 0.0, -1.0, 'T', 'N');  
      
      //-(G^T*R^-1*G)^-1*G^T*R^-1*dR_dthetai*R^-1*(Y-G*betaHat)
      //temp is a vector so solve only costs O(ntrend^2) ops
      solve_after_Chol_fact(dNbetaHat_dthetaN, Gtran_Rinv_G_Chol,temp);

      //temp3 is deps_dthetai
      matrix_mult(temp3, G, dNbetaHat_dthetaN, 0.0, -1.0);
      if(4<=maxObjDerMode) {
	//only store this if we're going to calculate the objective 
	//function's hessian
	deps_dtheta.putCols(temp3,i);
      }
      
      destVarianceMLE_dtheta(i)=(2.0*dot_product(rhs,temp3)-
                  		 dot_product(rhs,temp2))/numPoints;

      //this is correct _BECAUSE_ R^-1 and/or (and) dR_dthetai is/are _SYMMETRIC_
      trace_Rinv_dR_dthetai=dot_product(Rinv,dR_dthetai);

      matrix_mult(Gtran_Rinv_dR_thetai,Rinv_G,dR_dthetai,0.0,1.0,'T','N');
      matrix_mult(Gtran_Rinv_dR_thetai_Rinv_G,Gtran_Rinv_dR_thetai,Rinv_G,0.0,1.0,'N','N');
      trace_Gtran_Rinv_G_inv_Gtran_Rinv_dR_dthetai_Rinv_G=
	dot_product(Gtran_Rinv_G_inv,Gtran_Rinv_dR_thetai_Rinv_G);


      //We want the gradient of the negative of the log likelihood
      //gradObj(i)=0.5*(numPoints*destVarianceMLE_dtheta/estVarianceMLE+
      //	        trace_Rinv_dR_dthetai); 
      gradObj(i)=0.5*(destVarianceMLE_dtheta(i)/estVarianceMLE+
		      trace_Rinv_dR_dthetai/numPoints); //on a per point basis
    }

    prevObjDerMode=2;  //increase prevObjDerMode to the current value
    if((obj_der_mode==2)&&(con_der_mode<=prevConDerMode)) {
      //early exit because we have everything we need
      return;
    }
  }
  

  //******************************************************************
  //******************************************************************
  // calculating the hessian of the objective function starts here
  //******************************************************************
  //******************************************************************
  if((prevObjDerMode<4)&&(4<=obj_der_mode)) {
    hessObj.newSize(numTheta,numTheta);
    dR_dthetai.newSize(numPoints,numPoints);

    if(numTheta>1) {
      //consider this to be a temp variable, we do not need to retain its state
      //between sequential calls to this function only needed for j!=i
      dR_dthetaj.newSize(numPoints,numPoints);
    }

    dNbetaHat_dthetaN.newSize(ntrend); //This is another temporary variable 
    //(we do not need to retain its state between sequential calls to this 
    //function) but I gave it a name reflecting its contents. FOR 
    //CALCULATION OF THE HESSIAN OF THE OBJECTIVE FUNCTION it holds the 
    //mixed second derivative of the least squares coefficients with respect 
    //to theta(i) and theta(j).  FOR CALCULATION OF THE GRADIENT OF THE 
    //OBJECTIVE FUNCTION dNbetaHat_dthetaN holds the first derivative of the
    //vector of least squares coefficients with respect to theta(i)

    //consider this to be a temp variable, we do not need to retain its state
    //between sequential calls to this function
    d2eps_dthetai_dthetaj.newSize(numPoints);

    //consider this to be a temp variable, we do not need to retain its state
    //between sequential calls to this function
    d2R_dthetai_dthetaj.newSize(numPoints,numPoints);

    //consider this to be a temp variable, we do not need to retain its state
    //between sequential calls to this function
    dR_dthetai_Rinv_eps_minus_deps_dthetai.newSize(numPoints,numTheta);

    //consider this to be a temp variable, we do not need to retain its state
    //between sequential calls to this function
    Rinv_dR_dthetai_Rinv_eps_minus_deps_dthetai.newSize(numPoints,numTheta);

    temp.newSize(ntrend); //FOR THE OBJECTIVE'S HESSIAN this holds
    //G^T*R^-1*(-d2R/dtidtj*Rinv*eps
    //          +dR/dtj*R^-1*(dR/dti*R^-1*eps-deps/dti)
    //          +dR/dti*R^-1*(dR/dtj*R^-1*eps-deps/dtj))
    //it is something else for the objective's gradient. temp should always
    //be a vector of ntrend elements to avoid allocation and deallocation costs

    temp2.newSize(numPoints); //FOR THE OBJECTIVE'S HESSIAN this always holds
    //the j column of Rinv_dR_dthetai_Rinv_eps_minus_deps_dthetai, if we had 
    //data slicing we wouldn't have to copy this out, it is something else for
    //the objective's gradient. temp2 should already be the correct size

    temp3.newSize(numPoints);//FOR THE OBJECTIVE'S HESSIAN this holds
    //the j column of Rinv_dR_dthetai_Rinv_eps_minus_deps_dthetai, if we had 
    //data slicing we wouldn't have to copy this out, it is something else for
    //the objective's gradient

    temp4.newSize(numPoints);//always holds the j column of 
    //dR_dthetai_Rinv_eps_minus_deps_dthetai, if had data slicing we wouldn't
    //have to copy this out

    temp5.newSize(numPoints); //an honest to goodness temporary variable 
    //starts as temp5=-d2R/dtidtj*Rinv*eps and progresses to
    //temp5=-d2R/dtidtj*Rinv*eps
    //      +dR/dtj*R^-1*(dR/dti*R^-1*eps-deps/dti)
    //      +dR/dti*R^-1*(dR/dtj*R^-1*eps-deps/dtj)

    //a single loop over numPoints*numTheta elements should be slightly 
    //faster than nested loops over numTheta and numPoints
    int nk=numPoints*numTheta;
    for(int k=0; k<nk; ++k) {
      dR_dthetai_Rinv_eps_minus_deps_dthetai(k)=
	dR_dtheta_Rinv_eps(k)-deps_dtheta(k); //stored these when computed 
      //the gradient of the objective function
    }
    //precompute and store O(numPoints^2*numTheta) ops
    solve_after_Chol_fact(Rinv_dR_dthetai_Rinv_eps_minus_deps_dthetai,RChol,
			  dR_dthetai_Rinv_eps_minus_deps_dthetai);
    //solve_after_LU_fact(Rinv_dR_dthetai_Rinv_eps_minus_deps_dthetai,RLU,
    //		ipvt_RLU,dR_dthetai_Rinv_eps_minus_deps_dthetai);

    double d2est_variance_mle_dthetai_dthetaj;
    double inv_est_variance_mle=1.0/estVarianceMLE;
    
    //******************************************************************
    //by doing j in reverse order, these quantities for next i can be 
    //reused from previous j, but since we haven't looped over j yet we
    //have to set them for the first i (i=0);
    //do the diagonal separate, a few simplifcations to save work
    j=i=0; 
    dcorrMtx_dthetak(dR_dthetai, R, i);
    apply_nugget_eval(dR_dthetai); //apply nugget EVAL is correct, the additive term is constant with respect to theta so drops out
    Rinv_dR_dthetai_Rinv_eps_minus_deps_dthetai.getCols(temp2,i);
    dR_dthetai_Rinv_eps_minus_deps_dthetai.getCols(temp4,j);

    //the reasons I'm multiplying Rinv*dR_dthetai instead of doing 
    //LUsolve(R\dR_dthetai) are
    //1) we are only going to use this as part of a trace (in hessObj) so 
    //   there is less danger from round off error in the scalar quantity
    //2) it's significantly faster
    matrix_mult(Rinv_dR_dthetai,Rinv,dR_dthetai,0.0,1.0);
    //******************************************************************
    for(i=0; i<numTheta; ++i) {
      j=i; //do the diagonal separate, why 2 simplifcations to save work plus 
      //by doing j in reverse order can reuse some quanties from previous j 
      //(which was previous i+1, current i = previous i+1)

      //dcorrMtx_dthetak(dR_dthetaj, R, j); //don't need because of i=j simplification
      //apply_nugget_eval(dR_dthetai); //apply nugget EVAL is correct, the additive term is constant with respect to theta so drops out

      dcorrMtx_dthetak(d2R_dthetai_dthetaj, dR_dthetai, j); //nugget applied to dR_dthetai carries over to d2R_dthetai_dthetaj
      //Rinv_dR_dthetai_Rinv_eps_minus_deps_dthetai.getCols(temp3,j); //for i!=j

      //set temp5=-d2R/dtidtj*Rinv*eps
      matrix_mult(temp5,d2R_dthetai_dthetaj,rhs,0.0,-1.0);

      //set d2est_variance_mle_dthetai_dthetaj=
      //                         -eps^T*R^-1*d2R/dtidtj*R^-1*eps 
      //will have additional terms added below
      d2est_variance_mle_dthetai_dthetaj=dot_product(rhs,temp5);

      //before this line temp5=-d2R/dtidtj*Rinv*eps
      //set temp5=-d2R/dtidtj*Rinv*eps
      //          +dR/dtj*R^-1*(dR/dti*R^-1*eps-deps/dti)
      //          +dR/dti*R^-1*(dR/dtj*R^-1*eps-deps/dtj)
      //matrix_mult(temp5,dR_dthetai,temp3,1.0,1.0); //for i!=j
      //matrix_mult(temp5,dR_dthetaj,temp2,1.0,1.0); //for i!=j
      matrix_mult(temp5,dR_dthetai,temp2,1.0,2.0); //i=j simplification 

      //set temp=G^T*R^-1*(-d2R/dtidtj*Rinv*eps
      //                   +dR/dtj*R^-1*(dR/dti*R^-1*eps-deps/dti)
      //                   +dR/dti*R^-1*(dR/dtj*R^-1*eps-deps/dtj))
      matrix_mult(temp,Rinv_G,temp5,0.0,1.0,'T','N');

      //set d2betaHat/dtidtj=(G^T*R^-1*G)^-1*G^T*R^-1* 
      //                     (-d2R/dtidtj*Rinv*eps
      //                      +dR/dtj*R^-1*(dR/dti*R^-1*eps-deps/dti)
      //                      +dR/dti*R^-1*(dR/dtj*R^-1*eps-deps/dtj))
      //O(numPoints^2) ops
      solve_after_Chol_fact(dNbetaHat_dthetaN,Gtran_Rinv_G_Chol,temp);

      //set d2eps/dtidtj=-G*d2betaHat/dtidtj
      matrix_mult(d2eps_dthetai_dthetaj,G,dNbetaHat_dthetaN,0.0,-1.0);

      //before this line ...
      //d2estVarianceMLE_dthetai_dthetaj(i,j)=-eps^T*R^-1*d2R/dtidtj*R^-1*eps 
      //now set d2estVarianceMLE_dthetai_dthetaj(i,j) equal to
      //=(-eps^T*R^-1*d2R/dtidtj*R^-1*eps 
      //  +2.0*d2eps/dtidtj*R^-1*eps
      //  +2.0*(dR/dtj*p^-1*eps-deps/dtj)^T*R^-1*(dR/dti*p^-1*eps-deps/dti)
      // )/numPoints
      d2est_variance_mle_dthetai_dthetaj=
	(d2est_variance_mle_dthetai_dthetaj
	 +2.0*(dot_product(d2eps_dthetai_dthetaj,rhs)+
	       dot_product(temp4,temp2))
	 )/numPoints;

      //the sum of the dot products is the trace I was talking about above
      hessObj(i,j)=0.5*
	(-destVarianceMLE_dtheta(i)*inv_est_variance_mle
	 *destVarianceMLE_dtheta(j)*inv_est_variance_mle
	 +d2est_variance_mle_dthetai_dthetaj*inv_est_variance_mle
	 +(-dot_product(Rinv_dR_dthetai,Rinv_dR_dthetai) //i=j simplification
	   +dot_product(Rinv,d2R_dthetai_dthetaj)
	   )/numPoints); //on a per point basis


      //******************************************************************
      //do reverse order j so can reuse lowest j work for next i
      //******************************************************************
      for(j=numTheta-1; i<j; --j) {
	dcorrMtx_dthetak(dR_dthetaj, R, j);
	apply_nugget_eval(dR_dthetai); //apply nugget EVAL is correct, the additive term is constant with respect to theta so drops out

	dcorrMtx_dthetak(d2R_dthetai_dthetaj, dR_dthetai, j); //nugget applied to dR_dthetai carries over to d2R_dthetai_dthetaj
	Rinv_dR_dthetai_Rinv_eps_minus_deps_dthetai.getCols(temp3,j); //for i!=j

	//the reasons I'm multiplying Rinv*dR_dthetaj instead of doing 
	//LUsolve(R\dR_dthetaj) are
	//1) we are only going to use this as part of a trace (in hessObj) so 
	//   there is less danger from round off error in the scalar quantity
	//2) it's significantly faster
	matrix_mult(Rinv_dR_dthetaj,Rinv,dR_dthetaj,0.0,1.0);
	
	//set temp5=-d2R/dtidtj*Rinv*eps
	matrix_mult(temp5,d2R_dthetai_dthetaj,rhs,0.0,-1.0);

	//set d2est_variance_mle_dthetai_dthetaj=
	//                         -eps^T*R^-1*d2R/dtidtj*R^-1*eps 
        //will have additional terms added below
	d2est_variance_mle_dthetai_dthetaj=dot_product(rhs,temp5);

	//before this line temp5=-d2R/dtidtj*Rinv*eps
	//set temp5=-d2R/dtidtj*Rinv*eps
	//          +dR/dtj*R^-1*(dR/dti*R^-1*eps-deps/dti)
	//          +dR/dti*R^-1*(dR/dtj*R^-1*eps-deps/dtj)
	matrix_mult(temp5,dR_dthetai,temp3,1.0,1.0); //for i!=j
	matrix_mult(temp5,dR_dthetaj,temp2,1.0,1.0); //for i!=j

	//set temp=G^T*R^-1*(-d2R/dtidtj*Rinv*eps
	//                   +dR/dtj*R^-1*(dR/dti*R^-1*eps-deps/dti)
	//                   +dR/dti*R^-1*(dR/dtj*R^-1*eps-deps/dtj))
	matrix_mult(temp,Rinv_G,temp5,0.0,1.0,'T','N');

	//set d2betaHat/dtidtj=(G^T*R^-1*G)^-1*G^T*R^-1* 
	//                     (-d2R/dtidtj*Rinv*eps
	//                      +dR/dtj*R^-1*(dR/dti*R^-1*eps-deps/dti)
	//                      +dR/dti*R^-1*(dR/dtj*R^-1*eps-deps/dtj))
	//O(numPoints^2) ops
	solve_after_Chol_fact(dNbetaHat_dthetaN, Gtran_Rinv_G_Chol,temp);

	//set d2eps/dtidtj=-G*d2betaHat/dtidtj
	matrix_mult(d2eps_dthetai_dthetaj,G,dNbetaHat_dthetaN,0.0,-1.0);
	dR_dthetai_Rinv_eps_minus_deps_dthetai.getCols(temp4,j);
	//before this line ...
	//d2estVarianceMLE_dthetai_dthetaj=-eps^T*R^-1*d2R/dtidtj*R^-1*eps
	//now set d2estVarianceMLE_dthetai_dthetaj(i,j) equal to
	//=(-eps^T*R^-1*d2R/dtidtj*R^-1*eps 
	//  +2.0*d2eps/dtidtj*R^-1*eps
	//  +2.0*(dR/dtj*p^-1*eps-deps/dtj)^T*R^-1*(dR/dti*p^-1*eps-deps/dti)
	// )/numPoints
	d2est_variance_mle_dthetai_dthetaj=
	  (d2est_variance_mle_dthetai_dthetaj 
	   +2.0*(dot_product(d2eps_dthetai_dthetaj,rhs)+
		 dot_product(temp4,temp2))
	   )/numPoints;

	//the sum of the dot products is the trace I was talking about above
	hessObj(i,j)=0.5*
	  (-destVarianceMLE_dtheta(i)*inv_est_variance_mle
	   *destVarianceMLE_dtheta(j)*inv_est_variance_mle
	   +d2est_variance_mle_dthetai_dthetaj*inv_est_variance_mle
	   +(-dot_product(Rinv_dR_dthetai,Rinv_dR_dthetaj) //for i!=j 
	     +dot_product(Rinv,d2R_dthetai_dthetaj)
	     )/numPoints); //on a per point basis
	hessObj(j,i)=hessObj(i,j);
      } //for(j=numTheta-1; i<j; --j)

      if(i<numTheta-1) {
	//by doing j in reverse order, these quantities for next i (that is 
	//current i+1) can be reused from previous j (which is current i+1) 
	j=i+1;
	dR_dthetai.copy(dR_dthetaj); //nugget was already applied
	temp2.copy(temp3);
	Rinv_dR_dthetai_Rinv_eps_minus_deps_dthetai.getCols(temp2,i);
	Rinv_dR_dthetai.copy(Rinv_dR_dthetaj);
	//want next i=j temp4 to equal current i+1=j temp4 which, and since 
	//next i=currenty i+1 it already has the right value
      } //if(i<numTheta-1)
    } //for(i=0; i<numTheta; ++i)
   
    prevObjDerMode=2; 
    if((obj_der_mode==4)&&(con_der_mode<=prevConDerMode)) {
      //"early" exit because we have everything we need 
      //note that this should be the last loop but have this here to 
      //accomodate future expansion
      return;
    }
  }//if((prevObjDerMode<4)&&(4<=obj_der_mode))

  return;
}


void KrigingModel::getRandGuess(MtxDbl& guess) const
{
  int mymod = 1048576; //2^20 instead of 10^6 to be kind to the computer
  guess.newSize(1,numVarsr);
  double corr_length,tempdouble;
  for(int j=0;j<numVarsr;j++) {
    //tempdouble=(std::rand() % mymod);
    //tempdouble/=mymod;
    //printf("tempdouble=%g ",tempdouble);
    //corr_length=pow(2.0,tempdouble*(log2(max_corr_length)-log2(min_corr_length)) + log2(min_corr_length));    
      
    //corr_length=exp((std::rand() % mymod)*(maxNatLogCorrLen-minNatLogCorrLen)/mymod+
    //minNatLogCorrLen);
    //guess(j) = 1.0/(2.0*corr_length*corr_length); 
    guess(j) = (std::rand() % mymod)*(maxNatLogCorrLen-minNatLogCorrLen)/mymod+
      minNatLogCorrLen; //this returns a random nat_log_corr_len which is the space we need to search in

    //guess(j)=100.0;
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
  int chol_info;
  MtxDbl theta(1,numTheta);
  for(k=0; k<numTheta; ++k)
    theta(k)=0.5*exp(-2.0*nat_log_corr_len(k));

  R.newSize(numPoints,numPoints);
  
  temp.newSize(numPoints);
  correlation_matrix(theta); //assigns to member variable R

  if((ifChooseNug==true)||(nug<0.0))
    nug=0.0;

  apply_nugget_build();
  RChol.copy(R);
  double best_rcond;
  Chol_fact(RChol,chol_info,best_rcond);
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
      guess_theta(k)=0.5*exp(-2.0*guess(k));

    correlation_matrix(guess_theta); //assigns to member variable R
    apply_nugget_build();
    RChol.copy(R);
    Chol_fact(RChol,chol_info,rcond);
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

  if(best_rcond>rcond) {
    correlation_matrix(theta); //assigns to member variable R
    rcond=best_rcond;
  }

  MtxDbl best_theta(theta);

  //calculate the constraint that ensures that the correlation matrix is well conditioned
  //conceputally: condition_number < maxCondNum  
  //implemented as: largest_eigenvalue / maxCondNum - smallest_eigenvalue(s) < 0 (so we don't have to worry about dividing by zero when the smallest eigenvalue is zero and the scaling should be roughly consistant regardless of the number of points)
  MtxDbl all_eigvect(numPoints,numPoints), all_eigval(numPoints);

  eig_sym(all_eigvect, all_eigval, R);
  double small_eigval     = all_eigval(0); //smallest
  double large_eigval     = all_eigval(numPoints-1); //largest
  double constraint_value = (large_eigval+nug)/maxCondNum - (small_eigval+nug);

  if(constraint_value<=0){
    //printf("\n"); //for debug only
    return nat_log_corr_len;
  }
  
  MtxDbl small_eigvect(numPoints), large_eigvect(numPoints);
  MtxDbl dR_dthetak(numPoints,numPoints);
  MtxDbl grad_con(numTheta);
  double grad_con_dot_grad_con;


  double small_diff=0.00001;
  double prev_constraint_value=9999.9;
  int loopcount=0;
  int num_times_same_constraint=0;
  
  for(loopcount=0; loopcount<15; ++loopcount){
    
    if(fabs(constraint_value/prev_constraint_value-1.0)<small_diff)
      num_times_same_constraint++;
    else
      num_times_same_constraint=0;

    if((constraint_value<=0.0)||(num_times_same_constraint>=4))
	break;

    prev_constraint_value=constraint_value;
    //iteratively improve the condition number of the correlation matrix by taking a step in the opposite direction of the gradient of the constraint function with a step size that WOULD take the constraint to zero IF the constraint function was a linear function
    all_eigvect.getCols(small_eigvect,0);
    all_eigvect.getCols(large_eigvect,numPoints-1);
    grad_con_dot_grad_con=0.0;
    double dlarge_eigval_dthetak,dsmall_eigval_dthetak;
    double dlarge_eigval_dnat_log_corr_lenk, dsmall_eigval_dnat_log_corr_lenk;
    for(k=0; k<numTheta; k++) {
      dcorrMtx_dthetak(dR_dthetak, R, k);
    
      //caclulate the gradient of the constraint... 
      //using analytical derivatives of eigenvalues, the formula for 
      //exact arbitray order derivatives is found in 
      //  Jankovic, Exact nth Derivatives of Eigenvalues and Eigenvectors.
      //     Journal of Guidance, Control, and Dynamics, Vol 17. No 1,
      //     January-February 1994 pp 136-144
      //I (KRD) found this paper difficult to read but highly useful 
      //once it is deciphered.  To help with the deciphering, Eqns {1, 
      //20, and 22-24} are the only ones you will typically need for 
      //_linear_ algebra (example 1 is useful too).
      dlarge_eigval_dthetak=
	dot_product(large_eigvect,matrix_mult(temp,dR_dthetak,large_eigvect,0.0,1.0));
      dsmall_eigval_dthetak=
	dot_product(small_eigvect,matrix_mult(temp,dR_dthetak,small_eigvect,0.0,1.0));
      
      //grad_con(k)=dlarge_eigval_dthetak/maxCondNum-dsmall_eigval_dthetak;
      
      dlarge_eigval_dnat_log_corr_lenk=dlarge_eigval_dthetak*-2.0*theta(k);
      dsmall_eigval_dnat_log_corr_lenk=dsmall_eigval_dthetak*-2.0*theta(k);
      grad_con(k)=dlarge_eigval_dnat_log_corr_lenk/maxCondNum-dsmall_eigval_dnat_log_corr_lenk;

      if((nat_log_corr_len(k)<=minNatLogCorrLen)&&(grad_con(k)>0))
	grad_con(k)=0.0;
  
      grad_con_dot_grad_con+=grad_con(k)*grad_con(k);
    }
    
    

    double fraction=1.0, temp_double;
    for(k=0;k<numTheta;k++) {
      temp_double=(nat_log_corr_len(k)-minNatLogCorrLen)/(constraint_value/grad_con_dot_grad_con*grad_con(k));
      if(temp_double<0.0) {
	printf("temp_double=%g\n",temp_double);
	assert(!(temp_double<0.0));
      }
      if(temp_double<fraction)
	fraction=temp_double;
    }

    for(k=0;k<numTheta;k++) {
      //theta(k)-=constraint_value/grad_con_dot_grad_con*grad_con(k);
      nat_log_corr_len(k)-=fraction*constraint_value/grad_con_dot_grad_con*grad_con(k);
      theta(k)=0.5*exp(-2.0*nat_log_corr_len(k));
    }

    correlation_matrix(theta); //fills member variable R
    eig_sym(all_eigvect, all_eigval, R);
    small_eigval     = all_eigval(0);
    large_eigval     = all_eigval(numPoints-1);
    constraint_value = (large_eigval+nug)/maxCondNum - (small_eigval+nug);
  }

  if((constraint_value>0.0)&&(ifChooseNug==true)){
    nug=(large_eigval-maxCondNum*small_eigval)/(maxCondNum-1.0);
    //nug=constraint_value*maxCondNum/(maxCondNum-1.0); //would also work but is less clear
    if(nug>maxChooseNug)
      nug=maxChooseNug;      
  }

  //printf("makeGuessFeasible loopcount=%d constraint_value=%g nug=%g\n",loopcount,constraint_value,nug);
  
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
  if(constraintType.compare("eig")==0) //eigenvalue constraints can do analytical 1st derivatives
    assert(((maxObjDerMode==1)||(maxObjDerMode==3))&&((maxConDerMode==1)||(maxConDerMode==3))&&(maxConDerMode<=maxObjDerMode));
  else if(constraintType.compare("rcond")==0) //rcond constraint can NOT do analytical 1st derivatives
    assert(((maxObjDerMode==1)||(maxObjDerMode==3))&&(maxConDerMode==1));
  else
    assert((constraintType.compare("eig")==0)||(constraintType.compare("rcond")==0));

  if((maxObjDerMode==1)&&(maxConDerMode==1))
    opt.conminData.nfdg = 0; //use numerical  gradients of objective and constraints
  else if((maxObjDerMode==3)&&(maxConDerMode==3))
    opt.conminData.nfdg = 1; //use analytical gradients of objective and constraints 
  else if((maxObjDerMode==3)&&(maxConDerMode==1))
    opt.conminData.nfdg = 2; //uses analytical derivatives for the objective and numerical derivatives for constraints
  else
    assert(0);

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
