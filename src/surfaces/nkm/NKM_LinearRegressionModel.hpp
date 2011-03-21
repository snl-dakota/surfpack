#ifndef __LIN_REGRESS_MODEL_HPP__
#define __LIN_REGRESS_MODEL_HPP__

#include "NKM_SurfPack.hpp"
#include "NKM_SurfData.hpp"
#include "NKM_Optimize.hpp"
#include "NKM_SurfPackModel.hpp"

#include <cstdio>
#include <cstdlib>

namespace nkm {

/** 
    LinearRegressionModel: a class for creating and evaluating
    polynomial response surface models or (arbitrary?) order
*/
class LinearRegressionModel: public SurfPackModel
{

public:

  /// default constructor
  LinearRegressionModel()
  { /* empty constructor */ }
  
  /// Standard constructor
  LinearRegressionModel(const SurfData& sd);

  /* TODO: likely remove
  LinearRegressionModel(const LinearRegressionModel& other) : sd(other.sd) {*this=other; return;};
  
  LinearRegressionModel(const LinearRegressionProblem& problem,
			MtxDbl& EulAng_in);

  LinearRegressionModel(const LinearRegressionProblem& problem);



  inline LinearRegressionModel& operator=(const LinearRegressionModel& other){
    //sd=other.sd; 
    XR=other.XR; Y=other.Y; EulAng=other.EulAng; Rot=other.Rot; betaHat=other.betaHat; rms=other.rms; Poly=other.Poly; return *this;};

  */

  MtxDbl& evaluate_d1y(MtxDbl& d1y, const MtxDbl& xr) const {
    std::cerr << "evaluation of gradients not yet implemented in nkm::LinearRegressionModel" << std::endl;
    return d1y;
  }

  MtxDbl& evaluate_d2y(MtxDbl& d2y, const MtxDbl& xr) const {
    std::cerr << "evaluation of hessians not yet implemented in nkm::LinearRegressionModel" << std::endl;
    return d2y;
  }


  /// this logic is common to a significant number of "evalBasis" and
  /// "evaluate" functions so do it once and call it to make code
  /// maintenance easy
  MtxDbl& rot_xr(MtxDbl& xr_rot, const MtxDbl& rot_or_eul_ang, const MtxDbl& xr) const {
    return (rotate_xr(xr_rot,rot_or_eul_ang,xr));
  };

  /*{
    int nvarsr=xr.getNCols();
    if((rot_or_eul_ang.getNRows()==nvarsr)&&
       (rot_or_eul_ang.getNCols()==nvarsr)) {
      //rot_or_eul_ang is a rotation matrix
      matrix_mult(xr_rot,xr,rot_or_eul_ang,0.0,1.0);
    }
    else if((rot_or_eul_ang.getNRows()==((nvarsr*(nvarsr-1))/2))&&
	    (rot_or_eul_ang.getNCols()==1)) {
      //rot_or_eul_ang is a vector containing the Euler Angles for the rotation matrix
      MtxDbl rot; 
      gen_rot_mat(rot,rot_or_eul_ang,nvarsr);
      matrix_mult(xr_rot,xr,rot,0.0,1.0);
    }
    else{
      printf("Error in rot_xr(MtxDbl& xr_rot,MtxDbl& rot_or_eul_ang,MtxDbl& xr): rot_or_eul_ang has the wrong size!!!\n");
      assert(0);
    }
    return xr_rot;
  }
  */

  /** evaluates every polynomial basis function in "poly" at every
      point in xr, and returns them in matrix g in a reasonably
      efficient way; it works for arbitrary number of dimensions and
      polynomial order */
  MtxDbl& evalBasis(MtxDbl& g, const MtxInt& poly, const MtxDbl& xr) const
  {
    evaluate_poly_basis(g,poly,xr);
    return g;
  }

  /** evaluates g from a rotated xr (xr is preserved), this function
      wraps evalBasis(MtxDbl& g, MtxInt& poly, MtxDbl& xr) */
  MtxDbl& evalBasis(MtxDbl& g, const MtxInt& poly, const MtxDbl& rot_or_eul_ang, const MtxDbl& xr) const
  {
    MtxDbl xr_rot;
    rot_xr(xr_rot, rot_or_eul_ang, xr);
    evalBasis(g,poly,xr_rot);
    return g;
  }

  /** evaluates an instantiated linear regression models basis
      functions at point(s) xr AFTER rotating xr (the orginal xr is
      preserved) */
  MtxDbl& evalBasis(MtxDbl& g, const MtxDbl& xr) const
  {
    MtxDbl xr_rot;
    matrix_mult(xr_rot,xr,Rot,0.0,1.0);
    evalBasis(g,Poly,xr_rot);
    return g;
  }

  /// evaluates an instantiated linear regression model at a single point xr
  double evaluate(const MtxDbl& xr) const
  {
    MtxDbl g;
    evalBasis(g,xr);
    return (dot_product(g,betaHat));
  }

  /// evaluates an instantiated linear regression model at multiple points xr
  MtxDbl& evaluate(MtxDbl&y, const MtxDbl& xr) const
  {
    MtxDbl g;
    evalBasis(g,xr);
    matrix_mult(y,g,betaHat,0.0,1.0);
    return y;
  }

  /** evaluates a polynomial fit at a single point AFTER rotating xr
      (the original xr is preserved), it wraps evaluate(MtxInt& poly,
      MtxDbl& xr, MtxDbl& beta) */
  double evaluate(const MtxInt& poly, const MtxDbl& rot_or_eul_ang, 
		  const MtxDbl& xr, const MtxDbl& beta) const
  {
    MtxDbl xr_rot;
    rot_xr(xr_rot, rot_or_eul_ang, xr);
    return evaluate(poly,xr,beta);    
  }

  /// evaluation of a polynomial fit at a single point xr without rotating xr 
  double evaluate(const MtxInt& poly, const MtxDbl& xr, const MtxDbl& beta) const
  {
    assert((xr.getNRows()==1)&&(poly.getNRows()==beta.getNRows()));
    MtxDbl g; 
    evalBasis(g,poly,xr);
    return (dot_product(g,beta));
  }

  /** evaluation of a polynomial fit at multiple points xr AFTER
      rotating xr (the original xr is preserved), this function wraps
      evaluate(MtxDbl& y, MtxInt& poly, MtxDbl& xr, MtxDbl& beta) */
  MtxDbl& evaluate(MtxDbl& y, MtxInt& poly, 
		   MtxDbl& rot_or_eul_ang, MtxDbl& xr, 
		   MtxDbl& beta) {
    MtxDbl xr_rot;
    rot_xr(xr_rot,rot_or_eul_ang,xr);
    evaluate(y,poly,xr_rot,beta);    
    return y;
  }

  /// evaluation of a polynomial fit at multiple points xr without rotating xr
  MtxDbl& evaluate(MtxDbl& y, MtxInt& poly, MtxDbl& xr, MtxDbl& beta)
  {
    MtxDbl g;
    evalBasis(g,poly,xr);
    matrix_mult(y,g,beta,0.0,1.0); //this waprs BLAS
    return y;
  }

  int getNPoly()
  { return Poly.getNRows(); }

  int getNEulAng()
  { return EulAng.getNElems(); }

  double getRMS()
  { return rms; }

  /// the objective function is the RMS error of the least squares fit
  virtual double objective(const MtxDbl& EulerAng);

  //TODO: debating about whether or not to put an objective and gradient to determine the EulerAngles, this could be complicated, not convinced it would work well either 

private:

  SurfData sd;

  MtxDbl XR;

  MtxDbl Y;

  MtxInt Poly;

  MtxDbl EulAng;
  
  MtxDbl Rot;

  MtxDbl G;

  MtxDbl betaHat;

  double rms;

  int numVarsr;

  int numPoints; 

};

} // end namespace nkm

#endif
