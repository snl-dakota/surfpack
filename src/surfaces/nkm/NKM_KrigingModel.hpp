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
    MtxDbl corr_len(numTheta,1);
    for(int i=0; i<numTheta; ++i)
      corr_len(i,0)=std::exp(nat_log_corr_len(i,0));
    correlations.newSize(numTheta,1);
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
    MtxDbl corr_len(numTheta,1);
    for(int i=0; i<numTheta; ++i)
      corr_len(i,0)=std::exp(nat_log_corr_len(i,0));
    correlations.newSize(numTheta,1);
    get_theta_from_corr_len(correlations,corr_len);
    con_out.newSize(numConFunc,1);
    //MtxDbl theta(1,numTheta);
    for(int i=0; i<numTheta; ++i)
      correlations(i,0)=0.5*std::exp(-2.0*nat_log_corr_len(i,0));
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

  /// evaluate the trend function g(xr), using class member Poly
  inline MtxDbl& eval_trend_fn(MtxDbl& g, const MtxDbl& xr) {
    return (evaluate_poly_basis(g, flyPoly, Poly, xr));
  }

  inline MtxDbl& eval_der_trend_fn(MtxDbl& dg, const MtxInt& der, 
				   const MtxDbl& xr) {
    return (evaluate_poly_der_basis(dg, flyPoly, derivBetaHat, Poly, der, xr));
    return dg;
  }


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
  
  /** the order of the derivatives this Kriging Model was built for
      buildDerOrder=0  means function values only (regular Kriging)
      buildDerOrder=1  means function values plus gradients (Gradient Enhanced 
                       Kriging)
      buildDerOrder>=2 is currently not allowed, a later developer could 
                       implement Hessian Enhanced Kriging but KRD did not 
		       do this  */
  short buildDerOrder; 
  
  /** number of derivatives used to construct the Kriging/GEK model
      the zeroth-order derivative, i.e. function value itself counts
      as a derivative, so for Kriging nDer=1, for GEK nder=1+numVarsr */
  int nDer;

  /** a "Poly" style matrix (see "Poly" below) that stores derivatives 
      orders used to construct the Kriging/GEK model.  Der is a 
      numVarsr by nDer matrix. For Kriging Der is a numVarsr by 1 matrix 
      of zeros.  For GEK Der is a numVarsr by 1+numVarsr matrix whose 
      first (index zero) column is all zeros and columns with index 1 
      through numVarsr hold the identiy matrix which is the mixed partial
      derivative order representation of the gradient. */
  MtxInt Der;

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

  /** used to determine the "small feasible region" that we need to search to
      find good correlation lengths for the chosen correlation function */
  double aveDistBetweenPts;

  /// the upper bound of the small feasible region of correlation lengths
  double maxNatLogCorrLen;

  /// the lower bound of the small feasible region of correlation lengths
  double minNatLogCorrLen;

  /** the chosen natural log of the correlation LENGTHS (NOT correlation 
      PARAMETERS) */
  MtxDbl natLogCorrLen;

  /** the vector of correlation parameters (NOT correlation LENGTHS), these
      are the values determinened by the maximum likelihood optimization, the
      temporary in process version is called "theta" */
  MtxDbl correlations;

  /// NUMber of Real input VARiableS => NUMRVARS => NUMVARSR => numVarsr
  int numVarsr;

  /// NUMber of THETA for Real input variables... this is the number of correlation parameters (for real input variables), for Kriging numTheta=numVarsr, for radial basis functions numTheta=1 (radial basis functions, intended to be a derived class, are independent of direction)
  int numTheta;

  /** what optimization method did the user want us to use to determine the 
      correlation lengths of the Gaussian Process error model */
  std::string optimizationMethod;

  /// number of starting locations for (possibly multistart) local optimization to try
  int numStarts;

  /// maximum number of sets of roughness parameters to try
  int maxTrials;
  int maxTrialsGlobal; //used if optimization_method = global_local
  int maxTrialsLocal; //used if optimization_method = global_local

  ///used if optimization_method = global_local
  int maxTrialsGlobal; 

  ///used if optimization_method = global_local
  int maxTrialsLocal; 

  //"rcond" is now the only allowed constraint type, the eig option was removed
  //std::string constraintType;

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

  /** the trend basis functions evaluated at all available build data points
      in their original order.  It has npoly rows. For Kriging it has numPoints
      columns.  For GEK it has (1+numVarsr)*numPoints columns with each "whole
      point" appearing as as 1+numVarsr sequential columns */
  MtxDbl Gall;

  /** the transpose of the matrix of trend function evaluations at the 
      (likely reorderd) subset of points used to build the Kriging Model.  
      If GEK is used, it will also contain derivatives of the original 
      polynomial basis functions. The convention is that capital matrices 
      are for the data the model is built from, lower case matrices are 
      for arbitrary points to evaluate the model at */
  MtxDbl Gtran;

  /** true if the user said to use only the main effects (no interaction/mixed 
      terms) in the polynomial basis */
  bool ifReducedPoly; 

  /** this is what the user asked for, highest total order of any term in the 
      polynomial basis */
  int polyOrderRequested; 
  
  /** maxAllowedPolyOrder will equal polyOrderRequested unless the user
      didn't give us enough data, then it can be less than what they
      requested, will remove this once we start using Pivoted Cholesky
      to adaptively select an optimal subset of basis functions */
  int maxAllowedPolyOrder;

  /** this is what was actually used, can be less than what the user asked for 
      if the correlation matrix was ill-conditioned and we had to drop points 
      to fix the ill-conditioning and had to drop trend order to make it less 
      than the remaining points. */
  int polyOrder; 

  /** the number of equations needed for trend functions of order 0 through 
      polyOrderRequested */
  MtxInt numTrend; 

  /// the number of terms in the trend function that was actually used 
  int nTrend; 

  /** the polynomial powers for individual dimensions in each (trend) "basis 
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

  /// modified coefficients for on the fly derivatives of the trend functions
  MtxDbl derivBetaHat;

  /** the Z matrix, Z=Z(XR), 
      * for the Gaussian Correlation function
        Z(ij,k)=-(XR(i,k)-XR(j,k))^2
      * for the Exponential, Matern 1.5 and Matern 2.5 Correlation functions
        Z(ij,k)=-|XR(i,k)-XR(j,k)|
      * for the powered Exponential Correlation function (other than the 
        Exponential and Gaussian correlation functions)
        Z(ij,k)=-(|XR(i,k)-XR(j,k)|^powExpCorrFuncPow)
      The Z matrix facilitates the efficient evaluation of the correlation 
      matrix R and its derivatives with respect to theta (the vector of 
      correlation parameters);
      i.e. R=coefficient*exp(Z*theta) the coefficient is 1.0 for the powered
      Exponential family of correlation functions, it is not 1.0 for the 
      Matern 1.5 and Matern 2.5 correlation functions.  The convention is 
      that capital matrices are for the data the model is built from, 
      lower case matrices are for arbitrary points to evaluate model at 
      size(Z)=[numVarsr nchoosek(numPoints,2)] */
  MtxDbl Z; 

  /** working memory (so we don't constantly need to allocate and deallocate it)
      used during the calculation of R, equals Z^T*theta (matrix multiplication 
      is used), 
      size(Ztran_theta)=[nchoosek(numPoints,2) 1] */
  MtxDbl Ztran_theta;

  /** deltaXR(ij,k)=XR(i,k)-XR(j,k) for the strictly lower triangular part of R
      this is useful for efficient computation of derivative enhanced R matrix 
      for Gradient Enhanced Kriging, if Kriging is used this will be an
      empty matrix (i.e. space will not be allocated for it) 
      size(deltaXR)=[nchoosek(numPoints,2) numVarsr] 
      having a transposed order of Z should make it faster to compute the GEK R
      matrix sine it is "blocked" into (1+numVarsr) by (1+numVarsr) submatrices.
      The size of each submatrix is numPoints by numPoints (i.e. the size of 
      the Kriging R matrix) */
  MtxDbl deltaXR;

  /** the "correlation matrix," for either regular Kriging or Gradient Enhanced
      Kriging, after possible inclusion of a nugget, use of a nugget causes 
      the KrigingModel to smooth i.e. approximate (which is useful if you would
      like to account for known measurement noise) rather than interpolate and 
      can be used to fix ill-conditioning. Technically R is only a Correlation
      Matrix if all of it's diagonal elements = 1.0, this is not the case if
      you've added a nugget or if you're using Gradient Enhanced Kriging */
  MtxDbl R;

  //we will need Rinv if we want to evaluate the INTEGRAL of the adjusted variance; if that gets implemented we should compute Rinv as the "last step" of the construction of the Kriging/GEK model
  //MtxDbl Rinv; 
  //Rinv_Gtran*inv(G_Rinv_Gtran)*Rinv_Gtran^T would also be needed for analtyical integration of the adjusted variance, we should calculate this ONCE after the optimization in complete (when we clear stuff, or maybe we should just go ahead and calculate the integral (mean) of adjusted variance and store the answer in case anyone ever asks for it, then we wouldn't have to store Rinv or the longer matrix forever)

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

  /// part of infrastructure to allow masterObjectivesAndConstraints to just "return" (have copied out) the answer if the same point is used in sequential calls
  int prevObjDerMode;

  /// part of infrastructure to allow masterObjectivesAndConstraints to just "return" (have copied out) the answer if the same point is used in sequential calls
  int prevConDerMode;

  /// part of infrastructure to allow masterObjectivesAndConstraints to just "return" (have copied out) the answer if the same point is used in sequential calls
  MtxDbl prevTheta; //(numTheta,1)

  /// part of infrastructure to allow masterObjectivesAndConstraints to just "return" (have copied out) the answer if the same point is used in sequential calls
  int maxObjDerMode;

  /// part of infrastructure to allow masterObjectivesAndConstraints to just "return" (have copied out) the answer if the same point is used in sequential calls
  int maxConDerMode;

  /// the objective function for the optimization of correlation lengths, it's the negative "per equation" log likelihood function
  double obj;

  /// the vector (a numConFunc by 1 matrix) of constraint functions for the optimization of the correlation lengths, it only needs to be a vector for compatibility with the nkm::OptimizationProblem class, otherwise it could be a double
  MtxDbl con; 
};

} // end namespace nkm

#ifdef SURFPACK_HAVE_BOOST_SERIALIZATION
template< class Archive >
void nkm::KrigingModel::serialize(Archive & archive, 
			  const unsigned int version)
{  

  archive & boost::serialization::base_object<nkm::SurfPackModel>(*this);
  archive & buildDerOrder;
  archive & nDer;
  archive & Der;
  archive & corrFunc;
  archive & powExpCorrFuncPow;
  archive & maternCorrFuncNu;
  archive & aveDistBetweenPts;
  archive & maxNatLogCorrLen;
  archive & minNatLogCorrLen;
  archive & natLogCorrLen;
  archive & correlations;
  archive & numVarsr;
  archive & numTheta;
  archive & optimizationMethod;
  archive & ifUserSpecifiedCorrLengths;
  archive & numStarts;      
  archive & maxTrials;
  archive & maxTrialsGlobal;
  archive & maxTrialsLocal;
  archive & numConFunc;
  archive & maxCondNum;
  archive & ifChooseNug;
  archive & ifPrescribedNug;
  archive & nug;
  archive & iPtsKeep;
  archive & numPoints;
  archive & numPointsKeep;
  archive & numWholePointsKeep;
  archive & numExtraDerKeep;
  archive & numEqnAvail;
  archive & numRowsR;
  archive & ifHaveAnchorPoint;
  archive & iAnchorPoint;
  //don't archive XR since it's only a reference into SurfData .xr
  archive & XRreorder;
  //archive & Yall; //need this during the construction of a model but not afterward so don't archive
  archive & Y;
  //don't archive Gall we need this during the construction of a model but not afterward
  //don't archive Gtran we need this during the construction of a model but not afterward
  archive & ifReducedPoly;
  archive & polyOrderRequested;
  archive & maxAllowedPolyOrder; //will remove this once we've added the adaptive selection of a subset of basis functions
  archive & polyOrder;
  archive & numTrend;
  archive & nTrend;
  archive & Poly;
  //don't archive flyPoly, it's work space needed to efficiently evaluate a polynomial basis "on the fly" (it's a member variable only to avoid constant allocation and deallocation)
  archive & betaHat;
  //don't archive derivBeta, it's work space needed to efficiently evaluate derivatives of a polynomial (it's a member variable only to avoid constant allocation and deallocation)
  //don't archive Z, we need it during the construction of a model but not afterward
  //don't archive Ztran_theta, we need it during the construction of a model but not afterward
  //don't archive deltaXR, we need it during the construction of a model but not afterward
  //don't archive R, we need it during the construction of a model but not afterward, once we have ranking of candidate points by pivoteted cholesky we will use it as temporary variable space after the model is constructed but it still won't be something we want to retain
  //archive & Rinv; //not used, would be used for analytic integral of adjusted variance
  archive & RChol;
  //don't archive scaleRChol, we need it during the construction of a model but not afterward
  //don't archive sumAbsColR, we need it during the construction of a model but not afterward
  //don't archive oneNormR, we need it during the construction of a model but not afterward
  //don't archive lapackRcondR, weneed it during the construction of a model but not afterward
  //don't archive rcondDblWork, we need it during the construction of a model but not afterward
  //don't archive rcondIntWork, we need it during the construction of a model but not afterward
  archive & rcondR;
  archive & rcond_Gtran_Rinv_G;
  archive & Rinv_Gtran;
  archive & G_Rinv_Gtran_Chol;
  //don't archive G_Rinv_Gtran_Chol_Scale, we need it during the construction of a model but not afterward
  //don't archive G_Rinv_Gtran_Chol_DblWork, we need it during the construction of a model but not afterward
  //don't archive G_Rinv_Gtran_Chol_IntWork, we need it during the construction of a model but not afterward
  //don't archive G_Rinv_Y, we need it during the construction of a model but not afterward
  //don't archive eps, we need it during the construction of a model but not afterward
  archive & rhs;
  archive & estVarianceMLE;
  archive & likelihood;
  //don't archive prevObjDerMode, we need it during the construction of a model but not afterward
  //don't archive prevConDerMode, we need it during the construction of a model but not afterward
  //don't archive prevTheta, we need it during the construction of a model but not afterward
  archive & maxObjDerMode;
  archive & maxConDerMode;
  archive & obj;
  //don't archive con, we need it during the construction of a model but not afterward
}
BOOST_CLASS_EXPORT(nkm::KrigingModel)
#endif

#endif
