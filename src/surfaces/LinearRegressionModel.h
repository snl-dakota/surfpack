/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#ifndef __LINEAR_REGRESSION_MODEL_H__
#define __LINEAR_REGRESSION_MODEL_H__
#ifdef HAVE_CONFIG_H
#include "surfpack_config.h"
#endif

#include "surfpack_system_headers.h"
#include "SurfpackModel.h"
#include "SurfpackMatrix.h"
class SurfPoint;

class LRMBasisSet
{
public:
  VecVecUns bases;
  double eval(unsigned index, const VecDbl& x) const;
  double deriv(unsigned index, const VecDbl& x, const VecUns& vars) const;
  std::string asString() const;
  void add(const std::string& s_basis);
  unsigned size() const { return bases.size();}
};

class LinearRegressionModel : public SurfpackModel
{
public:
  LinearRegressionModel(const unsigned dims, const LRMBasisSet& bs_in, 
    const VecDbl& coeffs_in);
  virtual VecDbl gradient(const VecDbl& x) const;
  virtual std::string asString() const;
protected:
  virtual double evaluate(const VecDbl& x) const;
  LRMBasisSet bs;
  VecDbl coeffs;
friend class LinearRegressionModelTest;
};

struct Term {
  bool color;
  VecUns vars;
  Term(const VecUns& vars_in) : color(false), vars(vars_in) {}
};
  
///////////////////////////////////////////////////////////
///	Linear Regression Model Factory	
///////////////////////////////////////////////////////////

class LinearRegressionModelFactory : public SurfpackModelFactory 
{

public:
  LinearRegressionModelFactory();
  LinearRegressionModelFactory(const ParamMap& args);
  void setEqualityConstraints(unsigned asv,const SurfPoint& sp,  
    double valuePtr, VecDbl* gradientPtr, MtxDbl* hessianPtr);
  virtual SurfpackModel* Create(const SurfData& sd);
  virtual SurfpackModel* Create(const std::string& model_string);
  virtual void config();
  virtual unsigned minPointsRequired();
  VecDbl lrmSolve(const LRMBasisSet& bs, const ScaledSurfData& ssd);
  static LRMBasisSet CreateLRM(unsigned order, unsigned dims);
protected:
  unsigned order;
  MtxDbl eqConLHS;
  VecDbl eqConRHS;
};
#endif
