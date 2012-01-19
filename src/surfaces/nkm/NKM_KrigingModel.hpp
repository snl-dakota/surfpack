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

  std::string model_summary_string() const;

  // BMA TODO: can we redesign so these need not be public?
  void set_conmin_parameters(OptimizationProblem& opt) const;

  void set_direct_parameters(OptimizationProblem& opt) const;

  // Creating KrigingModels

  /// Default constructor
  KrigingModel() : ifChooseNug(false), nug(0.0), maxChooseNug(0.2), XR(sdBuild.xr), Y(sdBuild.y)
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
      correlations(0,i)=0.5*std::exp(-2.0*nat_log_corr_len(0,i));
    masterObjectiveAndConstraints(correlations, 1, 0);
    //printf("[objective]");
    return obj;
  };
  
  
  /* //useful because it shows how to convert the Hessian of the objective function from being with respect to theta to being with respect to natLogCorrLen, but the different (mountain) objective function didn't result a better model for the test cases and was considerably more expensive
  inline double objective(const MtxDbl& nat_log_corr_len) {
    //testing using "the volume of the mountain" (probability) rather than "the height of the mountain's peak" (likelihood=probability density) as the objective... it seems to work identically to maximum likelihood but require evaluation of higher derivatives so don't use should still print out hess_ob_out during model construction (ValidateMain tests) to verify that the hessian is working correctly 
    //the volume of an ellitpical paraboloid =(height/2)*prod(ellipse_radii)*volume of the unit hypersphere; each radius=1/(second derivative at endpoint of axis) (sign indicates direction), and those are the eigenvalues of the Hessian of the objective function. the 1/2 and volume of the hypersphere are constant with respect to theta, so they drop out.

    correlations.newSize(1,numTheta);
    //MtxDbl theta(1,numTheta);
    for(int i=0; i<numTheta; ++i)
      correlations(0,i)=0.5*std::exp(-2.0*nat_log_corr_len(0,i));
    masterObjectiveAndConstraints(correlations, 4, 0);

    //convert Hessian into log(correlation_length) space
    MtxDbl hess_obj_out(numTheta,numTheta);
    for(int i=0; i<numTheta; ++i)
      for(int j=0; j<numTheta; ++j)
	hess_obj_out(i,j)=hessObj(i,j)*4.0*correlations(0,i)*correlations(0,j);

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

    printf("currently you can't calculate analytical gradients of the objective function or constraints\n");
    assert(false);

    grad_obj_out.newSize(numTheta,1);
    correlations.newSize(1,numTheta);
    //MtxDbl theta(1,numTheta);
    for(int i=0; i<numTheta; ++i)
      correlations(0,i)=0.5*std::exp(-2.0*nat_log_corr_len(0,i));
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
    correlations.newSize(1,numTheta);
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
    for(int i=0; i<numTheta; ++i)
      correlations(0,i)=0.5*std::exp(-2.0*nat_log_corr_len(0,i));
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

  // BMA TODO: these docs need updating

  /** r(i,j)=exp(-sum_k theta(k) *(xr(i,k)-XR(j,k))^2); The convention is
      that capital matrices are for the data the model is built from, 
      lower case matrices are for arbitrary points to evaluate the model at */
  MtxDbl& correlation_matrix(MtxDbl& r, const MtxDbl& xr) const;

  /** dr(i,j)=d(r(i,j))/d(xr(i,k)) combining repeated calls can be used to 
      get arbitrary (mixed) higher order derivatives */
  MtxDbl& dcorrelation_matrix_dxI(MtxDbl& dr, const MtxDbl& r, 
				  const MtxDbl& xr, int Ider) const;
  MtxDbl& d2correlation_matrix_dxIdxK(MtxDbl& d2r, const MtxDbl& drI, const MtxDbl& r, const MtxDbl& xr, int Ider, int Kder) const;

  /* this function applies the nugget to an r matrix (not a member variable)
      in place (i.e. it changes r to p). i.e. it scales r by 1.0/(1.0+nug). 
      The convention is that capital matrices are for the data the model is 
      build from, lower case matrices are fore arbitrary points to evaluate 
      the model at */
  //MtxDbl& apply_nugget_eval(MtxDbl& r_to_p) const;

  /** R(i,j)=exp(-sum_k theta(k) *(XR(i,k)-XR(j,k))^2), where theta =
     corr_vec, implmented as R=exp(Z*theta) where Z=Z(XR),
     Z(ij,k)=-(XR(i,k)-XR(j,k))^2, */
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

  
  std::string optimizationMethod;

  /// number of starting locations for (possibly multistart) local optimization to try
  int numStarts;

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
  MtxDbl Rinv; //(numPoints,numPoints)

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

#endif
