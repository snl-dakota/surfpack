/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#ifndef __KRIGING_MODEL_HPP__
#define __KRIGING_MODEL_HPP__

//#include "surfpack_system_headers.h"
#include "NKM_SurfPack.hpp"
#include "NKM_SurfData.hpp"
#include "NKM_SurfPackModel.hpp"
#include "NKM_Optimize.hpp"
#include "NKM_LinearRegressionModel.hpp"
#include <map>
#include <string>

namespace nkm {

typedef std::map< std::string, std::string> ParamMap;

// enumerated type stored in nkm::KrigingModel::corrFunc see below for more details
enum {DEFAULT_CORR_FUNC, GAUSSIAN_CORR_FUNC, EXP_CORR_FUNC, POW_EXP_CORR_FUNC, MATERN_CORR_FUNC};

// BMA TODO: Use more descriptive names for variables?

/** 
    KrigingModel: a class for creating and evaluating Gaussian process
    emulators with constant, linear, or quadratic trend function.
    Options for:

    * choice of optimizer
    * coordinate rotation
    * evaluation with gradients
    * nugget to control ill-conditioning.
    * optimal subset selection to control ill-conditioning
    * correlation Function (powered exponential or matern families)
*/
class KrigingModel: public SurfPackModel
{

public:

  // MtxDbl& makeGuessFeasible(MtxDbl& nat_log_corr_len, OptimizationProblem *opt);

  std::string model_summary_string() const;

  // BMA TODO: can we redesign so these need not be public?
  void set_conmin_parameters(OptimizationProblem& opt) const;

  void set_direct_parameters(OptimizationProblem& opt) const;

  // Creating KrigingModels

  /// Default constructor
  KrigingModel() : ifChooseNug(false), nug(0.0), maxChooseNug(0.2), XR(sdBuild.xr) , Y(sdBuild.y)
  { /* empty constructor */ };
  
  /// Standard KrigingModel constructor
  KrigingModel(const SurfData& sd, const ParamMap& params);

  // BMA: in theory shouldn't need copy or assignment, given data
  // members below?!?

//   /// Copy constructor 
//   KrigingModel(const KrigingModel& other)
//   {  *this=other; assert(0);  } //effective C++ says not to use class assignment opperator in copy constructor

//   /// Assignment operator
//   KrigingModel& operator=(const KrigingModel& other)
//   { XR=other.XR; Y=other.Y; numVarsr=other.numVarsr; RLU=other.RLU; ipvt_RLU=other.ipvt_RLU; likelihood=other.likelihood; rhs=other.rhs; betaHat=other.betaHat; correlations=other.correlations; Poly=other.Poly; Rot=other.Rot; EulAng=other.EulAng; return *this; }

  /// After construction a Kriging model must be created with this
  /// function (TODO: add builtFlag for safety)
  void create();


  // Evaluating Kriging Models

  /// evaluate (y) the Kriging Model at a single point (xr is a Real row vector)
  virtual double evaluate(const MtxDbl& xr) const;

  /// evaluate (y) the Kriging Model at a collection of points xr, one per row
  virtual MtxDbl& evaluate(MtxDbl& y, const MtxDbl& xr) const;

  /// evaluate the KrigingModel's adjusted variance at a single point
  virtual double eval_variance(const MtxDbl& xr) const;

  /// evaluate the KrigingModel's adjusted variance at a collection of points xr, one per row
  virtual MtxDbl& eval_variance(MtxDbl& adj_var, const MtxDbl& xr) const;

  //double get_unadjusted_variance(){return (estVarianceMLE*scaler.unScaleFactorVarY());};
  
  /// evaluate the partial first derivatives with respect to xr of the models adjusted mean
  virtual MtxDbl& evaluate_d1y(MtxDbl& d1y, const MtxDbl& xr) const;

  /// evaluate the partial second derivatives with respect to xr of the models adjusted mean... this gives you the lower triangular, including diagonal, part of the Hessian(s), with each evaluation point being a row in both xr (input) and d2y(output)
  virtual MtxDbl& evaluate_d2y(MtxDbl& d2y, const MtxDbl& xr) const;

  // Helpers for solving correlation optimization problems

  /// the objective function, i.e. the negative log(likelihood);
  /// minimizing this produces a "KrigingModel" good)  
  inline double objective(const MtxDbl& nat_log_corr_len) {
    MtxDbl corr_len(1,numTheta);
    for(int i=0; i<numTheta; ++i)
      corr_len(0,i)=std::exp(nat_log_corr_len(0,i));
    correlations.newSize(1,numTheta);
    get_theta_from_corr_len(correlations,corr_len);
    masterObjectiveAndConstraints(correlations, 1, 0);
    //printf("[objective]");
    return obj;
  };
    
  /// the objective function, i.e. the negative log(likelihood), and
  /// its gradient; minimizing the objective function produces a good
  /// KrigingModel
  inline void objectiveAndGradient(double& obj_out, MtxDbl& grad_obj_out,
				   const MtxDbl& nat_log_corr_len) {

    printf("currently you can't calculate analytical gradients of the objective function or constraints\n");
    assert(false);

    grad_obj_out.newSize(numTheta,1);
    MtxDbl corr_len(1,numTheta);
    for(int i=0; i<numTheta; ++i)
      corr_len(0,i)=std::exp(nat_log_corr_len(0,i));
    correlations.newSize(1,numTheta);
    get_theta_from_corr_len(correlations,corr_len);
    masterObjectiveAndConstraints(correlations, 2, 0);
    obj_out=obj;
    //grad_obj_out.copy(gradObj);

    //convert from gradient with respect to theta to gradient with 
    //respect to nat_log_corr_len
    for(int i=0; i<numTheta; ++i)
      grad_obj_out(i,0)=gradObj(i,0)*-2.0*correlations(0,i);
    //printf("[grad_obj_out={%g",grad_obj_out(0,0));
    //for(int i=1; i<numTheta; ++i)
    //printf(", %g",grad_obj_out(i,0));
    //printf("}]");

    return;
  };
  
  /// objective plus condition number constraints
  //void objectiveAndConstraints(double& obj_out, MtxDbl& con_out, 
  inline void objectiveAndConstraints(double& obj_out, MtxDbl& con_out, 
				      const MtxDbl& nat_log_corr_len) {
    //printf("entered objectiveAndConstraints\n");  fflush(stdout);
    MtxDbl corr_len(1,numTheta);
    for(int i=0; i<numTheta; ++i)
      corr_len(0,i)=std::exp(nat_log_corr_len(0,i));
    correlations.newSize(1,numTheta);
    get_theta_from_corr_len(correlations,corr_len);
    con_out.newSize(numConFunc,1);
    //MtxDbl theta(1,numTheta);
    for(int i=0; i<numTheta; ++i)
      correlations(0,i)=0.5*std::exp(-2.0*nat_log_corr_len(0,i));
    //printf("about to enter masterObjectiveAndConstraints\n"); fflush(stdout);
    masterObjectiveAndConstraints(correlations, 1, 1);
    //printf("left masterObjectiveAndConstraints\n"); fflush(stdout);
    obj_out=obj;
    for(int i=0; i<numConFunc; i++){
      //printf("i=%d ",i); fflush(stdout);
      con_out(i,0)=con(i,0);
    }
    //con_out.copy(con);
    //printf("[objectiveAndConstraints]");
    //printf("leaving objectiveAndConstraints\n");  fflush(stdout);
    return;
  };

  /// objective plus condition number constraints with gradients
  inline void objectiveAndConstraintsAndGradients(double& obj_out, 
						  MtxDbl& con_out, 
						  MtxDbl& grad_obj_out, 
						  MtxDbl& grad_con_out, 
						  const MtxDbl& nat_log_corr_len) {
    printf("currently you can't calculate analytical gradients of the objective function or constraints\n");
    assert(false);

    con_out.newSize(numConFunc,1);
    grad_obj_out.newSize(numTheta,1);
    grad_con_out.newSize(numConFunc,numTheta);
    //MtxDbl theta(1,numTheta);
    MtxDbl corr_len(1,numTheta);
    for(int i=0; i<numTheta; ++i)
      corr_len(0,i)=std::exp(nat_log_corr_len(0,i));
    correlations.newSize(1,numTheta);
    get_theta_from_corr_len(correlations,corr_len);
    masterObjectiveAndConstraints(correlations, 2, 2);
    obj_out=obj;
    for(int i=0; i<numConFunc; ++i)
      con_out(i,0)=con(i,0);

    //convert from gradient with respect to theta to gradient with 
    //respect to nat_log_corr_len
    for(int j=0; j<numTheta; ++j) {
      grad_obj_out(j,0)=gradObj(j,0)*-2.0*correlations(0,j);
      for(int i=0; i<numConFunc; ++i)
	grad_con_out(i,j)=gradCon(i,j)*-2.0*correlations(0,j);
    }

    //printf("[grad_obj_out={%g",grad_obj_out(0));
    //for(int i=1; i<numTheta; ++i)
    //printf(", %g",grad_obj_out(i));
    //printf("}]");

    return;
  };

  /// return the Number of Trend functions, the trend is represented by an
  /// arbitrary order multidimensional polynomial, individual trend functions
  /// are the separate additive terms in that multidimensional polynomial
  inline int getNTrend() const
  { return (Poly.getNRows());   } 

  // return the likelihood of this model
  inline double getLikelihood()
  { return likelihood; }

  static int min_coefficients(int nvars, int poly_order) 
  {
    return num_multi_dim_poly_coef(nvars,poly_order)+nvars;
  };

  void getRandGuess(MtxDbl& guess) const;

private:
  
#ifdef SURFPACK_HAVE_BOOST_SERIALIZATION
  // allow serializers access to private data
  friend class boost::serialization::access;
  /// serializer for derived class SurfPoint data
  template<class Archive> 
  void serialize(Archive & archive, const unsigned int version);
#endif

  // helper functions
  void preAllocateMaxMemory();
  void equationSelectingCholR();
  void nuggetSelectingCholR();

  /// this function calculates the objective function (negative log
  /// likelihood) and/or the constraint functions and/or their analytical
  /// gradients and/or the hessian of the objective function using a 
  /// precompute and store (store across sequential calls to this function) 
  /// strategy to reduce the computational cost.  To ensure that precomputed
  /// and stored values are not changed externally this function return no
  /// output, instead member variables must be copied out by wrapper functions.
  /// The objective and contraint derivative modes are bit flags, i.e. 
  /// each is the sum of 2^(all orders of derivative you want). KRD 2010.05.13
  void masterObjectiveAndConstraints(const MtxDbl& theta, int obj_der_mode, 
				     int con_der_mode);

  //void set_conmin_parameters(OptimizationProblem& opt) const;

  /// evaluate the trend function g(xr), using specified {Poly, Rot}
  MtxDbl& eval_trend_fn(MtxDbl& g, const MtxInt& poly, 
			const MtxDbl& rot_or_eul_ang, const MtxDbl& xr) const;

  /// evaluate the trend function g(xr), using class members {Poly, Rot}
  MtxDbl& eval_trend_fn(MtxDbl& g, const MtxDbl& xr) const;

  //the following matern_1pt5_... and matern_2pt5_... functions don't really need to be member functions as they don't access any data members

  /// multiply exponential corr func by this to get matern 1.5 corr func
  inline double matern_1pt5_coef(double theta_abs_dx) const{
    return 1.0+theta_abs_dx;
  };
  /// multiply matern 1.5 corr func by this to get d1 of matern 1.5 corr func
  inline double matern_1pt5_d1_mult_r(double theta, double dx) const{
    return -theta*theta*dx/matern_1pt5_coef(theta*std::fabs(dx));
  };
  /** multiply matern 1.5 corr func by this to get d2 of matern 1.5 corr func
      1D MATERN_CORR_FUNC 1.5 r(x1,x2) is twice+ differential except 
      where x1==x2 this is correct for x1!=x2 */
  inline double matern_1pt5_d2_mult_r(double theta, double dx) const{
    return theta*theta*(1.0-2.0/matern_1pt5_coef(theta*std::fabs(dx)));
  };


  /// multiply exponential corr func by this to get matern 2.5 corr func
  inline double matern_2pt5_coef(double theta_abs_dx) const{
    return 1.0+theta_abs_dx+theta_abs_dx*theta_abs_dx/3.0;
  };
  /// multiply matern 2.5 corr func by this to get d1 of matern 2.5 corr func
  inline double matern_2pt5_d1_mult_r(double theta, double dx) const{
    double theta_abs_dx=theta*std::fabs(dx);
    return -theta*theta*dx*(1-theta_abs_dx)/(3*matern_2pt5_coef(theta_abs_dx));
  };
  /// multiply matern 2.5 corr func by this to get d2 of matern 2.5 corr func
  inline double matern_2pt5_d2_mult_r(double theta, double dx) const{
    double theta_abs_dx=theta*std::fabs(dx);
    return -theta*theta*
      (1-(2.0/3.0+2.0*theta_abs_dx)/matern_2pt5_coef(theta_abs_dx));
  };


  // BMA TODO: these docs need updating

  /** converts from correlation lengths to theta
      for powered exponential (including exponential and Gaussian)
          theta=1/(powExpCorrLenPow*corr_len^powExpCorrLenPow)
      for matern (excluding Gaussian)
          theta= sqrt(2*maternCorrFuncNu)/corr_len  */
  MtxDbl& get_theta_from_corr_len(MtxDbl& theta, const MtxDbl& corr_len) const;
  /** converts from theta to correlation lengths 
      for powered exponential (including exponential and Gaussian)
          theta=1/(powExpCorrLenPow*corr_len^powExpCorrLenPow)
      for matern (excluding Gaussian)
          theta= sqrt(2*maternCorrFuncNu)/corr_len  */
  MtxDbl& get_corr_len_from_theta(MtxDbl& corr_len, const MtxDbl& theta) const;

  /** r(i,j)=corr_func(xr(i,:),XR(j,:);theta(:)) choices for correlation 
      function are gaussian, exponential, powered exponential with 1<power<2, 
      and matern with nu=1.5 or 2.5 (gaussian and exponential are pulled out
      for efficient implementation and because they belong to both families).
      The convention is that capital matrices are for the data the model is 
      built from, lower case matrices are for arbitrary points to evaluate 
      the model at */
  MtxDbl& correlation_matrix(MtxDbl& r, const MtxDbl& xr) const;

  /** dr(i,j)=d(r(i,j))/d(xr(i,k)) combining repeated calls can be used to 
      get arbitrary (mixed) higher order derivatives */
  MtxDbl& dcorrelation_matrix_dxI(MtxDbl& dr, const MtxDbl& r, 
				  const MtxDbl& xr, int Ider) const;
  MtxDbl& d2correlation_matrix_dxIdxK(MtxDbl& d2r, const MtxDbl& drI, const MtxDbl& r, const MtxDbl& xr, int Ider, int Kder) const;

  /** R(i,j)=corr_func(XR(i,:),XR(j,:);theta(:)) where choices for for 
      correlation function are gaussian, exponential, powered exponential 
      with 1<power<2, and matern with nu=1.5 or 2.5 (gaussian and exponential 
      are pulled out for efficient implementation and because they belong 
      to both families).  All correlation functions are implemented as
      R=something.*exp(Z*theta) (with reshapes) the something depends on the
      correlation function (for gaussian, exponential, and powered exponential
      that something is 1) for the matern function that something is
      matern_1pt5_coef or matern_2pt5_coef) of course the definition of Z 
      and theta differs for different correlation functions.  Note that Z 
      only stores what is needed to compute the strictly lower (BELOW the 
      diagonal) part of R to save memory and computation.  R is symmetric 
      and has ones on the diagonal. */
  void correlation_matrix(const MtxDbl& corr_vec);

  /** this function applies the nugget to the R matrix (a member variable)
      and stores the result in R (another member variable), i.e. it adds 
      nug to the diagonal of R. The convention is that capital matrices 
      are for the data the model is built from, lower case matrices are for 
      arbitrary points to evaluate the model at.  Once the nugget is added
      R is not strictly a correlation matrix */
  void apply_nugget_build();

  /** the Z matrix, Z=Z(XR), its definitition depends on the correlation 
      function
          for the gaussian correlation function 
              Z(ij,k)=-(XR(i,k)-XR(j,k))^2, 
          for the exponetial and matern correlation functions
              Z(ij,k)=-|XR(i,k)-XR(j,k)|     
	  for the powered exponential correlation function
	      Z(ii,k)=-|XR(i,k)-XR(j,k)|^powExpCorrFuncPow
	      where 1<powExpCorrFuncPow<2
      the Z matrix facilitates the efficient evaluation of the correlation 
      matrix R... R=something.*exp(Z*theta) 
      note that Z only holds what is needed to compute the strictly lower
      (below the diagonal) portion of R to save memory and computation.
      R is symmetric and has ones on the diagonal.  The convention is that 
      capital matrices are for the data the model is built from, lower 
      case matrices are for arbitrary points to evaluate the model at, 
      the Z and XR matrices are member variables so they don't need to be 
      passed in */
  MtxDbl& gen_Z_matrix();
  
  /** stores which correlation function we using, major choices are 
          powered exponential with 1<=power<=2 and 
          matern with nu=0.5,1.5,2.5 or "infinity" 
      There are 2 special cases that belong to both families and are 
      pulled out for efficient implementation these are
          gaussian correlation function
              equals powered exponential with power=2 
              equals matern with nu="infinity"
          exponential correlation function 
              equals powered exponential with power=1
	      equals matern with nu = 0.5
      after these are pulled out we have 
      powered exponential with 1<power<2 and
      matern with nu = 1.5 or 2.5
      corrFunc stores an enumerated type
  */
  short corrFunc;
  double powExpCorrFuncPow;
  double maternCorrFuncNu;
  
  std::string optimizationMethod;

  /// number of starting locations for (possibly multistart) local optimization to try
  int numStarts;

  /// maximum number of sets of roughness parameters to try
  int maxTrials;
  int maxTrialsGlobal; //used if optimization_method = global_local
  int maxTrialsLocal; //used if optimization_method = global_local

  /// should contain either "eig" or "rcond"
  std::string constraintType;

  /** ifChooseNug=1 tells KrigingModel to choose the smallest nugget it 
      needs to fix ill conditioning, ifChooseNug=0 tells KrigingModel not 
      to choose one (the default value for the nugget is zero) but the 
      user still has the option to prescribe a nugget the Nugget must be 
      a positive number */
  bool ifChooseNug; 

  /// which if any preset nugget formula should be used to calculate a nugget size.
  int nuggetFormula;

  bool ifUserSpecifiedCorrLengths;
  
  int numVarsr;

  /// NUMber of THETA for Real input variables... this is the number of correlation parameters (for real input variables), for Kriging numTheta=numVarsr, for radial basis functions numTheta=1 (radial basis functions, a derived class, are independent of direction)
  int numTheta;

  double aveDistBetweenPts;
  double maxNatLogCorrLen;
  double minNatLogCorrLen;

  int numPoints; 
  int numRowsR;
  
  /** the number of constraint FUNTIONS (typically these are nonlinear), 
      this number does NOT include box edge constraints for the inputs,
      those are handled separately */
  int numConFunc; 

  /** the nugget value sets the ammount of smoothing (approximation instead
      of interpolation) that the KrigingModel will use, it can also be used
      to fix ill conditioning, setting ifChooseNug=1 tells KrigingModel to 
      choose the smallest nugget needed to fix ill conditioning */
  double nug;

  /** the maximum value the KrigingModel is allowed to choose to fix ill
      conditioning, however the user can prescribe a larger value to nug.
      if nug=maxChooseNug is not sufficient to decrease the condtion number
      of the correlation matrix below maxCondNum, then the correlation 
      parameter will be chosen to reduce R's condition number to the 
      prescribed level */
  double maxChooseNug;

  /** the correlation matrix R is considered to be "ill-conditioned" if
      if it's condition number exceeds this value.  We are using the term
      "ill-conditioned" loosely, because the qualtity of a KrigingModel
      is frequently improved when the condition number is restricted to be 
      less than the number of points (yes that is a heuristic and we're 
      looking for a better one) */
  double maxCondNum;

  /** the input the model was constructed from; convention is capital
      matrices are data model is built from, lower case matrices are
      arbitrary points to evaluate model at, using XR instead of X in
      anticipation of mixed real and integer input
  */
  MtxDbl& XR;

  /** the output the model was constructed from; convention is capital
      matrices are data model is built from, lower case matrices are
      arbitrary points to evaluate model at */
  MtxDbl Y;   

  /** the polynomial powers for individual dimensions in each "basis 
      function" (for now the only choice of polynomial basis functions are
      multidimensional monomials) in a multidimensional polynomial of 
      arbitrary order */
  MtxInt Poly;  
  bool ifReducedPoly; /// only use main effects (no interaction/mixed terms) in the polynomial basis
  int polyOrderRequested; /// this is what the user asked for, highest total order of any term in the polynomial basis
  int polyOrder; /// this is what was actually used, can be less than what the user asked for if the correlation matrix was ill conditioned and we had to drop points to fix the ill conditioning and had to drop trend order to make it less than the remaining points.
  MtxInt numTrend; //the number of equations needed for trend functions of order 0 through polyOrderRequested
  int nTrend;

  /** the Euler angles to generate the input dimensions' Rotation matrix */
  MtxDbl EulAng; 

  MtxDbl Rot; ///the input dimensions' Rotation matrix

  /** The lower triangular part of the Cholesky decomposition of R (the 
      correlation matrix after possible modification by the inclusion of 
      a nugget).  Keep this around to evaluate the adjusted variance. 
      The convention is that capital matrices are for the data the model 
      is built from, lower case matrices are for arbitrary points to 
      evaluate the model at */
  MtxDbl RChol;
  MtxDbl scaleRChol;
  MtxDbl sumAbsColR;
  MtxDbl oneNormR;
  MtxDbl lapackRcondR;
  MtxDbl rcondDblWork;
  MtxInt rcondIntWork;
  MtxInt iEqnKeep;
  bool ifHaveAnchorPoint;
  int  iAnchorPoint;
  int numEqnAvail;
  int numEqnKeep;

  MtxDbl Yall;
  MtxDbl Gall;

  /** LU decomposition of R (the correlation matrix after possible
      modification by the inclusion of a nugget).  Keep this around to 
      evaluate the adjusted variance. We may want to replace with Rinv, 
      that would facilitate efficient evaluation of the integral of 
      adjusted variance;  The convention is that capital matrices are for 
      the data the model is built from, lower case matrices are for 
      arbitrary points to evaluate the model at */
  //MtxDbl RLU; 

  /** pivot indices for LU decomposition of R (the possibly modified by
      a nugget correlation matrix R).  Keep this around to
      evaluate the adjusted variance.  We may want to replace with Rinv,
      that would facilitate efficient evaluation of the integral of
      adjusted variance */
  //MtxInt ipvt_RLU; 

  /// the log(likelihood) of the Kriging Model
  double likelihood; 

  /** rhs=Rinv*(Y-G(XR)*betaHat); The convention is that capital matrices 
      are for the data the model is built from, lower case matrices are 
      for arbitrary points to evaluate the model at */
  MtxDbl rhs; 

  /// the vector of coefficients of the trend functions (unadjusted mean)
  MtxDbl betaHat;

  /// the vector of correlation parameters (NOT correlation LENGTHS)
  MtxDbl correlations;

  /// the natural log of the corrletion LENGTHS (NOT correlation PARAMETERS)
  MtxDbl natLogCorrLen;

  /** the rcond (estimated reciprocal of the condition number) of the 
      modified correlation matrix, R */
  double rcondR; 
  double rcond_Gtran_Rinv_G;

  /** the matrix of trend function evaluations, G=G(XR); The convention is
      that capital matrices are for the data the model is built from, lower
      case matrices are for arbitrary points to evaluate the model at */
  MtxDbl G;

  /** the Z matrix, Z=Z(XR), Z(ij,k)=-(XR(i,k)-XR(j,k))^2, facilitates
      the efficient evaluation of the correlation matrix R and its
      derivatives with respect to theta (the vector of correlation
      parameters); The convention is that capital matrices are for the
      data the model is built from, lower case matrices are for 
      arbitrary points to evaluate model at */
  MtxDbl Z;
  MtxDbl Ztheta;

  /// misc members moved from KrigingProblem?!?

  /** the correlation matrix, after possible inclusion of a
      nugget, use of a nugget causes the KrigingModel to smooth i.e.
      approximate rather than interpolate and can be used to fix 
      ill-conditioning. */
  MtxDbl R;

  ///keep around to evaluate the integral of adjusted variance, this currently only getting calculated when optimization_method=local is selected
  //MtxDbl Rinv; //(numPoints,numPoints)

  /// keep around to evaluate adjusted variance
  MtxDbl Rinv_G;

  /// keep around to evaluate adjusted variance
  //MtxDbl Gtran_Rinv_G_LU;
  MtxDbl Gtran_Rinv_G_Chol;
  MtxDbl Gtran_Rinv_G_Chol_Scale;
  MtxDbl Gtran_Rinv_G_Chol_DblWork;
  MtxInt Gtran_Rinv_G_Chol_IntWork;

  /// keep around to evaluate adjusted variance
  double estVarianceMLE;

  MtxDbl temp;
  MtxDbl temp2;


  int prevObjDerMode;
  int prevConDerMode;
  MtxDbl prevTheta; //(1,numTheta)


  int maxObjDerMode;
  int maxConDerMode;
  double obj;
  MtxDbl gradObj; //(numTheta);
  MtxDbl hessObj; //(numTheta,numTheta);
  MtxDbl con; //(numConFunc);
  MtxDbl gradCon; //(numConFunc,numTheta);
};

} // end namespace nkm

#ifdef SURFPACK_HAVE_BOOST_SERIALIZATION
template< class Archive >
void nkm::KrigingModel::serialize(Archive & archive, 
			  const unsigned int version)
{  

  archive & boost::serialization::base_object<nkm::SurfPackModel>(*this);
  archive & corrFunc;
  archive & powExpCorrFuncPow;
  archive & maternCorrFuncNu;
  archive & optimizationMethod;          
  archive & numStarts;      
  archive & maxTrials;
  archive & maxTrialsGlobal;
  archive & maxTrialsLocal;
  archive & constraintType;
  archive & ifChooseNug;
  archive & nuggetFormula;
  archive & ifUserSpecifiedCorrLengths;
  archive & numVarsr;
  archive & numTheta;
  archive & aveDistBetweenPts;
  archive & maxNatLogCorrLen;
  archive & minNatLogCorrLen;
  archive & numPoints;
  archive & numRowsR;
  archive & numConFunc;
  archive & nug;
  archive & maxChooseNug;
  archive & maxCondNum;
  archive & XR;
  archive & Y;
  archive & Poly;
  archive & ifReducedPoly;
  archive & polyOrderRequested;
  archive & polyOrder;
  archive & numTrend;
  archive & nTrend;
  archive & EulAng; //need to keep Rot, but EulAng is more human readable, but I should probably scrap both of them and not use polynomials of rotated X as a trrend function (the infra structure is there but it's not being used and I'll need to remove it if I am to calculate the integral of the adjusted variance)
  archive & Rot; //see previous comment
  archive & RChol;
  archive & scaleRChol;
  //archive & sumAbsColR; //need this during the constuction of a model but not afterward
  //archive & oneNormR; //need this during the constuction of a model but not afterward
  //archive & lapackRcondR; //need this during the constuction of a model but not afterward
  //archive & rcondDblWork; //need this during the constuction of a model but not afterward
  //archive & rcondIntWork; //need this during the constuction of a model but not afterward
  archive & iEqnKeep;
  archive & ifHaveAnchorPoint;
  archive & iAnchorPoint;
  archive & numEqnAvail;
  archive & numEqnKeep;
  //archive & Yall; //need this during the constuction of a model but not afterward
  //archive & Gall; //need this during the constuction of a model but not afterward
  archive & likelihood;
  archive & rhs;
  archive & betaHat;
  archive & correlations;
  archive & natLogCorrLen;
  archive & rcondR;
  archive & rcond_Gtran_Rinv_G;
  //archive & G; //need this during the constuction of a model but not afterward
  //archive & Z; //need this during the constuction of a model but not afterward
  //archive & Ztheta; //need this during the constuction of a model but not afterward
  //archive & R; //need this during the constuction of a model but not afterward
  //archive & Rinv; //not used, would be used for integral of adjusted variance
  //Rinv_G*inv(Gtran_Rinv_G)*Rinv_G^T would also be needed should calculate this ONCE after the optimization in complete (when we clear stuff, or maybe we should just go ahead and calculate the integral (mean) of adjusted variance and store the answer in case anyone ever asks for it, then we wouldn't have to store Rinv or the longer matrix forever)
  archive & Rinv_G;
  archive & Gtran_Rinv_G_Chol;
  archive & Gtran_Rinv_G_Chol_Scale;
  //archive & Gtran_Rinv_G_Chol_DblWork; //need this during the constuction of a model but not afterward
  //archive & Gtran_Rinv_G_Chol_IntWork; //need this during the constuction of a model but not afterward
  archive & estVarianceMLE;
  //archive & temp; //need this during the constuction of a model but not afterward
  //archive & temp2; //need this during the constuction of a model but not afterward
  //archive & prevObjDerMode; //need this during the constuction of a model but not afterward
  //archive & prevConDerMode; //need this during the constuction of a model but not afterward
  //archive & prevTheta; //need this during the constuction of a model but not afterward
  archive & maxObjDerMode;
  archive & maxConDerMode;
  archive & obj;

  //archive & gradObj; //not needed because analytical derivatives removed
  //archive & hessObj; //not needed because analytical derivatives removed
  //archive & con; //not needed because analytical derivatives removed
  //archive & gradCon; //not needed because analytical derivatives removed
}
BOOST_CLASS_EXPORT(nkm::KrigingModel)
#endif

#endif
