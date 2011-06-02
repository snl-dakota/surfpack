#include <cmath>
#include <iostream>
#include "NKM_SurfPackModel.hpp"
#include "NKM_LinearRegressionModel.hpp"
#include "NKM_SurfPack.hpp"

namespace nkm {

using std::cout;
using std::endl;

//#define __LIN_REGRESS_ERR_CHECK__

/***********************************************************************/
/***********************************************************************/
/**** Unit Test functions for LinearRegression start here           ****/
/***********************************************************************/
/***********************************************************************/

/***********************************************************************/
/***********************************************************************/
/**** Unit Test functions for LinearRegression end here             ****/
/***********************************************************************/
/***********************************************************************/

/***********************************************************************/
/***********************************************************************/
/**** LinearRegressionModelFactory member functions start here      ****/
/***********************************************************************/
/***********************************************************************/

/***********************************************************************/
/***********************************************************************/
/**** LinearRegressionModelFactory member functions end here        ****/
/***********************************************************************/
/***********************************************************************/

/***********************************************************************/
/***********************************************************************/
/**** LinearRegressionModel member functions start here             ****/
/***********************************************************************/
/***********************************************************************/

/** constructor for when no vector of rotation angles is specified,
    i.e. don't rotate, i.e. rotation angles of zero, i.e. rotate by the
    identity matrix */
LinearRegressionModel::
LinearRegressionModel(const SurfData& sd_in)
  :sd(sd_in), numVarsr(sd.getNVarsr()), numPoints(sd.getNPts()), XR(sd.xr)
{
  sd.getY(Y);

  // TODO: allow user specification
  const bool estimate_euler_angles = true;

  // from problem
  assert(numPoints>0);

  // number of design variables, should optimization later be used
  int num_design_var = nchoosek(numVarsr,2); //nchoosek(numVarsr,2);
  //TODO: use experimentation to determine the number of guesses
  //printf("LRP: num_guess=%d\n",num_guess);
  // number of guesses should optimization later be used
  int num_guess = num_design_var*10; 

  int npoly;
  int poly_order_in = -99999;
  if(poly_order_in!=-99999) {
#ifdef __LIN_REGRESS_ERR_CHECK__
    assert(poly_order_in>=0);
#endif
    npoly = nchoosek(poly_order_in + numVarsr, poly_order_in);
    if(npoly>numPoints) {
      printf("Error in LinearRegressionProblem()!!! You have requested more basis functions than there are points (i.e. G^T*G will be singular)!!! Choose a lower order polynomial to fit. For least squares we recommend using 2*npoly<npts where npoly=nchoosek(poly_order+nvars,poly_order);\n");
      assert(npoly<=numPoints);
    }    
  }
  else{
    poly_order_in=0; 
    npoly=1;
    while(2*npoly<=numPoints){
      poly_order_in++;
      npoly=nchoosek(poly_order_in + numVarsr, poly_order_in);      
    }
    poly_order_in--;
    npoly=nchoosek(poly_order_in + numVarsr, poly_order_in);            
  }
  if(poly_order_in<=1) {
    //rotation only makes a difference if you're using second or higher order polynomials
    //printf("reset num_guess=1\n");
    num_guess = 1;
  }
  //printf("\nnumPoints=%d numVarsr=%d poly_order=%d npoly=%d\n",numPoint,numVarsr,poly_order_in,npoly);
  //printf("LRP numGuess=5%d\n",numGuess);
  Poly.newSize(npoly, numVarsr);
  multi_dim_poly_power(Poly, numVarsr, poly_order_in);

  if (estimate_euler_angles) {

    // setup data for optimization problem

    //problem.getInitGuess(eul_ang);
    
    OptimizationProblem opt(*this, num_design_var);
    
    // set bounds on angles to estimate
    double pi=2.0*acos(0.0);
    // upperbounds 0=2pi rest=pi puts axis one equally in all octants,
    // all =pi puts some axis equally in all octants
    for (int j=0; j<num_design_var; ++j) {
      opt.lower_bound(j, 0.0);
      opt.upper_bound(j, pi);
      opt.initial_iterate(j,0.0);
      // TODO: restore this capability?  Don't convey special meaning with zero
      //initGuess.zero();  ///all zeros means don't rotate the axis
    }

    opt.best_guess_optimize(num_guess);

    EulAng = opt.best_point();
    // initialize Rot
    gen_rot_mat(Rot,EulAng,XR.getNCols());

  }
  else {
    EulAng.newSize(1, nchoosek(numVarsr, 2)); 
    EulAng.zero();
    Rot.newSize(numVarsr, numVarsr);
    Rot.identity();
  }
  
  // code common to constructors (regardless of euler angle spec)
  evalBasis(G, XR);
  least_squares(G, betaHat, Y);

  MtxDbl eps = Y;
  matrix_mult(eps, G, betaHat, 1.0, -1.0); //eps_new=eps_old-G*betaHat;

  rms = eps(0,0)*eps(0,0);
  for(int i=1; i<numPoints; ++i)
    rms += eps(i,0)*eps(i,0);
  rms = std::sqrt(rms);

}

/// the objective function of a LinearRegressionProblem is the rms
/// error of the polynomial fit, it depends on the input coordinates
/// rotation matrix's Euler Angles
double LinearRegressionModel::objective(const MtxDbl& euler_angle) {
  //printf("numVarsr=%d",numVarsr); fflush(stdout);
  //printf(" euler_angle.size=[%d %d]",euler_angle.getNRows(),euler_angle.getNCols()); fflush(stdout);
  //printf(" eulerangle(0,0)=%12.6g\n",euler_angle(0,0));fflush(stdout);

  // TODO: consider local varible here now that classes merged
  gen_rot_mat(Rot, euler_angle, numVarsr);

  MtxDbl xr; 
  matrix_mult(xr,XR,Rot,0.0,1.0);
  LinearRegressionModel::evalBasis(G,Poly,xr);
  //do an unconstrained least squares, add constrainment in later
  least_squares(G, betaHat, Y);

  //calculate the error
  MtxDbl eps = Y;
  matrix_mult(eps, G, betaHat, 1.0, -1.0); //eps_new=eps_old-G*betaHat;
  rms = eps(0,0)*eps(0,0);
  for(int i=1; i<numPoints; ++i)
    rms += eps(i,0)*eps(i,0);
  rms = std::sqrt(rms);
  return rms;
}

/*
///evalBasis evaluates every polynomial basis function in "poly" at every point in xr, and returns them in matrix g in a reasonably efficient way; it works for arbitrary number of dimensions and polynomial order
MtxDbl& LinearRegressionModel::evalBasis(MtxDbl& g, MtxInt& poly, MtxDbl& xr) 
{
  int numPoints=xr.getNRows();
  int numVarsr=xr.getNCols();
  int npoly=poly.getNRows();
#ifdef __LIN_REGRESS_ERR_CHECK__
  if(!((numVarsr>0)&&(numVarsr==poly.getNCols())))
    assert((numVarsr>0)&&(numVarsr==poly.getNCols()));
#endif
  g.newSize(numPoints,npoly);
  
  int ivar, ipoly, ipt;
  for(ipoly=0; ipoly<npoly; ipoly++) {
#ifdef __LIN_REGRESS_ERR_CHECK__
    for(ivar=0;ivar<numVarsr;ivar++)
      assert(poly(ipoly,ivar)>=0);
#endif
    
    //save multiplications by not doing anything for polynomial powers of zero
    for(ivar=0;ivar<numVarsr;ivar++)
      if(poly(ipoly,ivar)!=0) 
	break;

    if(ivar==numVarsr) {
      //if all polynomial powers are zero then prod(anything^0)=1.0
      for(ipt=0; ipt<numPoints; ipt++)
	g(ipt,ipoly)=1.0;
    }
    else {

      //do simple assignment for the first non-zero power dimension (save assignment of 1.0 and the product of 1.0 times the power of xr)
      switch(poly(ipoly,ivar)) {
      case 0:
	//we shouldn't need this but better safe than sorry
	break;
      case 1:
	//linear functions are so simple and common that we don't want to use the "expensive" pow function
	for(ipt=0; ipt<numPoints; ipt++)
	  g(ipt,ipoly)=xr(ipt,ivar);
	break;
      case 2:
	//quadratic functions are so simple and common that we don't want to use the "expensive" pow function
	for(ipt=0; ipt<numPoints; ipt++)
	  g(ipt,ipoly)=xr(ipt,ivar)*xr(ipt,ivar);
	break;
      default:
	//we pulled out cubics as being simple and common too, but they're not as simple or common, typical use case for Kriging is power<=2, and we have to stop somewhere
	for(ipt=0; ipt<numPoints; ipt++)
	  g(ipt,ipoly)=std::pow(xr(ipt,ivar),poly(ipoly,ivar));
      }
      ivar++;
     
      //do every dimension after the first non-zero power dimension in the general fashion... i.e. we have to multiply these by the non-unity product we already have
      for(; ivar<numVarsr; ivar++) 
	switch(poly(ipoly,ivar)) {
	case 0:
	  //we DO need this
	  break;
	case 1:
	  //linear functions are so simple and common that we don't want to use the "expensive" pow function
	  for(ipt=0; ipt<numPoints; ipt++)
	    g(ipt,ipoly)*=xr(ipt,ivar);
	  break;
	case 2:
	  //quadratic functions are so simple and common that we don't want to use the "expensive" pow function
	  for(ipt=0; ipt<numPoints; ipt++)
	    g(ipt,ipoly)*=xr(ipt,ivar)*xr(ipt,ivar);
	  break;
	default:
	  //we pulled out cubics as being simple and common too, but they're not as simple or common, the typical use case for Kriging is power<=2, and we have to stop somewhere
	  for(ipt=0; ipt<numPoints; ipt++)
	    g(ipt,ipoly)*=std::pow(xr(ipt,ivar),poly(ipoly,ivar));
	}
    }
  }
      
  return g;
}
*/


/***********************************************************************/
/***********************************************************************/
/**** LinearRegressionModel member functions end here               ****/
/***********************************************************************/
/***********************************************************************/

} // end namespace nkm
