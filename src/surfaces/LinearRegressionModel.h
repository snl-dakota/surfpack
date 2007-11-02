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

//typedef std::vector<double> VecDbl;
//typedef std::vector<double>::const_iterator VecDblIt;
//typedef std::vector<unsigned> VecUns;
//typedef std::vector<unsigned>::const_iterator VecUnsIt;
//typedef std::vector< std::vector< unsigned > > VecVecUns;

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
//Factory type stuff
  static LinearRegressionModel Create(const SurfData& sd);
  static LRMBasisSet CreateLRM(unsigned order, unsigned dims);
  static VecDbl lrmSolve(const LRMBasisSet& bs, const ScaledSurfData& ssd);
// End factory stuff
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
  
#endif
