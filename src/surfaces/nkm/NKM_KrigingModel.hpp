/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#ifndef __KRIGING_MODEL_HPP__
#define __KRIGING_MODEL_HPP__
//#ifdef HAVE_CONFIG_H
//#include "surfpack_config.h"
//#endif

//#include "surfpack_system_headers.h"
//#include "SurfpackModel.h"
#include "NKM_SurfPack.hpp"
#include "NKM_SurfData.hpp"
#include "NKM_SurfPackModel.hpp"
#include "NKM_Optimize.hpp"
#include "NKM_LinearRegressionModel.hpp"
#include <map>
#include <string>

namespace nkm {

typedef std::map< std::string, std::string> ParamMap;

// BMA TODO: Use more descriptive names for variables?

/** 
    KrigingModel: a class for creating and evaluating Gaussian process
    emulators with constant, linear, or quadratic trend function.
    Options for:

    * choice of optimizer
    * coordinate rotation
    * evaluation with gradients
    * nugget to control ill-conditioning.
*/
class KrigingModel: public SurfPackModel
{

public:

  // BMA TODO: can we redesign so these need not be public?
  void set_conmin_parameters(OptimizationProblem& opt) const;

  void set_direct_parameters(OptimizationProblem& opt) const;

  // Creating KrigingModels

  /// Default constructor
  KrigingModel() : XR(sdBuild.xr), Y(sdBuild.y), ifChooseNug(false), nug(0.0), maxChooseNug(0.2)
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
  virtual MtxDbl& evaluate(MtxDbl& y, const MtxDbl& xr);

  /// evaluate the KrigingModel's adjusted variance at a single point
  double eval_variance(const MtxDbl& xr);

  /// evaluate the KrigingModel's adjusted variance at a collection of points xr, one per row
  MtxDbl& eval_variance(MtxDbl& adj_var, const MtxDbl& xr);

  MtxDbl& evaluate_d1y(MtxDbl& d1y, const MtxDbl& xr);

  MtxDbl& evaluate_d2y(MtxDbl& d2y, const MtxDbl& xr);

  // Helpers for solving correlation optimization problems

  /// adjust correlations to be feasible with respect to condition
  /// number constraints
  MtxDbl& makeGuessFeasible(MtxDbl& nat_log_corr_len, 
			    OptimizationProblem *opt);

  /// the objective function, i.e. the negative log(likelihood);
  /// minimizing this produces a "good" KrigingModel)
  
  
  inline double objective(const MtxDbl& nat_log_corr_len) {
    correlations.newSize(1,numTheta);
    //MtxDbl theta(1,numTheta);
    for(int i=0; i<numTheta; ++i)
      correlations(i)=0.5*exp(-2.0*nat_log_corr_len(i));
    masterObjectiveAndConstraints(correlations, 1, 0);
    //printf("[objective]");
    return obj;
  };
  
  
  /*
  inline double objective(const MtxDbl& nat_log_corr_len) {
    //testing using "the volume of the mountain" (probability) rather than "the height of the mountain's peak" (likelihood=probability density) as the objective... it seems to work identically to maximum likelihood but require evaluation of higher derivatives so don't use should still print out hess_ob_out during model construction (ValidateMain tests) to verify that the hessian is working correctly 
    //the volume of an ellitpical paraboloid =(height/2)*prod(ellipse_radii)*volume of the unit hypersphere; each radius=1/(second derivative at endpoint of axis) (sign indicates direction), and those are the eigenvalues of the Hessian of the objective function. the 1/2 and volume of the hypersphere are constant with respect to theta, so they drop out.

    correlations.newSize(1,numTheta);
    //MtxDbl theta(1,numTheta);
    for(int i=0; i<numTheta; ++i)
      correlations(i)=0.5*exp(-2.0*nat_log_corr_len(i));
    masterObjectiveAndConstraints(correlations, 4, 0);

    //convert Hessian into log(correlation_length) space
    MtxDbl hess_obj_out(numTheta,numTheta);
    for(int i=0; i<numTheta; ++i)
      for(int j=0; j<numTheta; ++j)
	hess_obj_out(i,j)=hessObj(i,j)*4.0*correlations(i)*correlations(j);

    MtxDbl eig_vals(numTheta); //all the eigenvalues of a Hessian should be positive
    eig_sym(eig_vals,hess_obj_out);
    double test_obj=obj;
    for(int i=0; i<numTheta; ++i)
      test_obj/=eig_vals(i);
	  
    return test_obj;
  };
  */


  /// the objective function, i.e. the negative log(likelihood), and
  /// its gradient; minimizing the objective function produces a good
  /// KrigingModel
  inline void objectiveAndGradient(double& obj_out, MtxDbl& grad_obj_out,
				   const MtxDbl& nat_log_corr_len) {

    grad_obj_out.newSize(numTheta);
    correlations.newSize(1,numTheta);
    //MtxDbl theta(1,numTheta);
    for(int i=0; i<numTheta; ++i)
      correlations(i)=0.5*exp(-2.0*nat_log_corr_len(i));
    masterObjectiveAndConstraints(correlations, 2, 0);
    obj_out=obj;
    //grad_obj_out.copy(gradObj);

    //convert from gradient with respect to theta to gradient with 
    //respect to nat_log_corr_len
    for(int i=0; i<numTheta; ++i)
      grad_obj_out(i)=gradObj(i)*-2.0*correlations(i);
    //printf("[grad_obj_out={%g",grad_obj_out(0));
    //for(int i=1; i<numTheta; ++i)
    //printf(", %g",grad_obj_out(i));
    //printf("}]");

    return;
  };
  
  /// objective plus condition number constraints
  //void objectiveAndConstraints(double& obj_out, MtxDbl& con_out, 
  inline void objectiveAndConstraints(double& obj_out, MtxDbl& con_out, 
				      const MtxDbl& nat_log_corr_len) {
    //printf("entered objectiveAndConstraints\n");  fflush(stdout);
    correlations.newSize(1,numTheta);
    con_out.newSize(numConFunc);
    //MtxDbl theta(1,numTheta);
    for(int i=0; i<numTheta; ++i)
      correlations(i)=0.5*exp(-2.0*nat_log_corr_len(i));
    //printf("about to enter masterObjectiveAndConstraints\n"); fflush(stdout);
    masterObjectiveAndConstraints(correlations, 1, 1);
    //printf("left masterObjectiveAndConstraints\n"); fflush(stdout);
    obj_out=obj;
    for(int i=0; i<numConFunc; i++){
      //printf("i=%d ",i); fflush(stdout);
      con_out(i)=con(i);
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
    con_out.newSize(numConFunc);
    grad_obj_out.newSize(numTheta);
    grad_con_out.newSize(numConFunc,numTheta);
    //MtxDbl theta(1,numTheta);
    for(int i=0; i<numTheta; ++i)
      correlations(i)=0.5*exp(-2.0*nat_log_corr_len(i));
    masterObjectiveAndConstraints(correlations, 2, 2);
    obj_out=obj;
    for(int i=0; i<numConFunc; ++i)
      con_out(i)=con(i);

    //convert from gradient with respect to theta to gradient with 
    //respect to nat_log_corr_len
    for(int j=0; j<numTheta; ++j) {
      grad_obj_out(j)=gradObj(j)*-2.0*correlations(j);
      for(int i=0; i<numConFunc; ++i)
	grad_con_out(i,j)=gradCon(i,j)*-2.0*correlations(j);
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
  { return (Poly.getNRows()); //XR.getNCols()+1); 
	    } 

  // return the likelihood of this model
  inline double getLikelihood()
  { return likelihood; }

  static int min_coefficients(int nvars, int poly_order) 
  {
    return num_multi_dim_poly_coef(nvars,poly_order)+nvars;
  };

  void getRandGuess(MtxDbl& guess) const;

private:

  // helper functions

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

  // BMA TODO: these docs need updating

  /** r(i,j)=exp(-sum_k theta(k) *(xr(i,k)-XR(j,k))^2); The convention is
      that capital matrices are for the data the model is built from, 
      lower case matrices are for arbitrary points to evaluate the model at */
  MtxDbl& correlation_matrix(MtxDbl& r, const MtxDbl& xr) const;

  /** dr(i,j)=d(r(i,j))/d(xr(i,k)) combining repeated calls can be used to 
      get arbitrary (mixed) higher order derivatives */
  MtxDbl& dcorrelation_matrix_dxk(MtxDbl& dr, const MtxDbl& r, 
				  const MtxDbl& xr, int k);

  /** this function applies the nugget to an r matrix (not a member variable)
      in place (i.e. it changes r to p). i.e. it scales r by 1.0/(1.0+nug). 
      The convention is that capital matrices are for the data the model is 
      build from, lower case matrices are fore arbitrary points to evaluate 
      the model at */
  MtxDbl& apply_nugget_eval(MtxDbl& r_to_p) const;

  /** R(i,j)=exp(-sum_k theta(k) *(XR(i,k)-XR(j,k))^2), where theta =
     corr_vec, implmented as R=exp(Z*theta) where Z=Z(XR),
     Z(ij,k)=-(XR(i,k)-XR(j,k))^2, */
  //MtxDbl 
  void correlation_matrix(const MtxDbl& corr_vec);

  /** this function applies the nugget to the R matrix (a member variable)
      and stores the result in R (another member variable), i.e. it adds 
      nug to the diagonal of R, scales the result by 1.0/(1.0+nug), and
      stores that result in R. The convention is that capital matrices are 
      for the data the model is built from, lower case matrices are for 
      arbitrary points to evaluate the model at */
  void apply_nugget_build();

  //static MtxDbl corrMtx(MtxDbl& R, int nrowsXR, MtxDbl& Z, MtxDbl& corr_vec);

  /** the Z matrix, Z=Z(XR), Z(ij,k)=-(XR(i,k)-XR(j,k))^2, facilitates
      the efficient evaluation of the correlation matrix
      R... R=exp(Z*theta)... and its derivatives with respect to theta
      (the vector of correlation
      parameters)... dR_dthetak=-Z(:,k).*R(:) (MATLAB notation);
      The convention is that capital matrices are for the data the model 
      is built from, lower case matrices are for arbitrary points to 
      evaluate the model at, the Z and XR matrices are member variables so 
      they don't need to be passed in */
  MtxDbl& gen_Z_matrix();

  /** dR_dthetak=-Z(:,k).*R(:) (MATLAB notation), and 
      d2R_dthetai_dthetaj=-Z(:,j).*dR_dthetai; The convention is that
      capital matrices are for the data the model is built from, lower
      case matrices are for arbitrary points to evaluate the model at. 
      Z is a member variables so it doesn't need to be passed in, R needs 
      to be passed in so we can evaluate second derivatives with the same
      function
  */
  MtxDbl& dcorrMtx_dthetak(MtxDbl& dR_dthetak, const MtxDbl& R, const int k);

  // data

  
  std::string optimizationMethod;

  /// maximum number of sets of roughness parameters to try
  int maxTrials;

  /// should contain either "eig" or "rcond"
  std::string constraintType;

  /** ifChooseNug=1 tells KrigingModel to choose the smallest nugget it 
      needs to fix ill conditioning, ifChooseNug=0 tells KrigingModel not 
      to choose one (the default value for the nugget is zero) but the 
      user still has the option to prescribe a nugget the Nugget must be 
      a positive number */
  bool ifChooseNug; 

  bool ifUserSpecifiedCorrLengths;
  
  int numVarsr;

  /// NUMber of THETA for Real input variables... this is the number of correlation parameters (for real input variables), for Kriging numTheta=numVarsr, for radial basis functions numTheta=1 (radial basis functions, a derived class, are independent of direction)
  int numTheta;


  double aveDistBetweenPts;
  double maxNatLogCorrLen;
  double minNatLogCorrLen;

  int numPoints; 
  
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
  MtxDbl& Y;   

  /** the polynomial powers for individual dimensions in each "basis 
      function" (for now the only choice of polynomial basis functions are
      multidimensional monomials) in a multidimensional polynomial of 
      arbitrary order */
  MtxInt Poly;
  
  /** the Euler angles to generate the input dimensions' Rotation matrix */
  MtxDbl EulAng; 

  MtxDbl Rot; //the input dimensions' Rotation matrix

  /** The lower triangular part of the Cholesky decomposition of R (the 
      correlation matrix after possible modification by the inclusion of 
      a nugget).  Keep this around to evaluate the adjusted variance. 
      The convention is that capital matrices are for the data the model 
      is built from, lower case matrices are for arbitrary points to 
      evaluate the model at */
  MtxDbl RChol;

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

  ///keep around to evaluate the integral of adjusted variance
  MtxDbl Rinv; //(numPoints,numPoints)

  /// keep around to evaluate adjusted variance
  MtxDbl Rinv_G;

  /// keep around to evaluate adjusted variance
  //MtxDbl Gtran_Rinv_G_LU;
  MtxDbl Gtran_Rinv_G_Chol;

  /// keep around to evaluate adjusted variance
  //MtxInt ipvt_Gtran_Rinv_G_LU;

  double estVarianceMLE;


  MtxDbl dNbetaHat_dthetaN; //(ntrend)
  MtxDbl temp; //(ntrend)
  MtxDbl temp2; //(numPoints)
  MtxDbl temp3; //(numPoints)
  MtxDbl temp4; //(numPoints)
  MtxDbl temp5; //(numPoints)
  MtxDbl d2eps_dthetai_dthetaj; //(numPoints)
  //MtxDbl P; //(numPoints,numPoints)
  MtxDbl dR_dthetai; //(numPoints,numPoints)
  MtxDbl dR_dthetaj; //(numPoints,numPoints)
  MtxDbl Rinv_dR_dthetai; //(numPoints,numPoints)
  MtxDbl Rinv_dR_dthetaj; //(numPoints,numPoints)  
  MtxDbl d2R_dthetai_dthetaj; //(numPoints,numPoints)
  MtxDbl dR_dthetai_Rinv_eps_minus_deps_dthetai; //(numPoints,numTheta)
  MtxDbl Rinv__dR_dthetai_Rinv_eps_minus_deps_dthetai; //(numPoints,numTheta)

  
  int prevObjDerMode;
  int prevConDerMode;
  MtxDbl prevTheta; //(1,numTheta)
  MtxDbl allEigVect; //(numPoints,numPoints)
  MtxDbl allEigVal; //(numPoints)
  //MtxDbl Z; //(numPoints*numPoints,numTheta)
  //MtxDbl G; //(numPoints,ntrend)
  MtxDbl deps_dtheta; //(numPoints,numTheta)
  MtxDbl dR_dtheta_Rinv_eps; //(numPoints,numTheta)
  MtxDbl destVarianceMLE_dtheta; //(numTheta)

  int maxObjDerMode;
  int maxConDerMode;
  double obj;
  MtxDbl gradObj; //(numTheta);
  MtxDbl hessObj; //(numTheta,numTheta);
  MtxDbl con; //(numConFunc);
  MtxDbl gradCon; //(numConFunc,numTheta);
};

} // end namespace nkm

#endif
