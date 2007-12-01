/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#ifndef __DIRECT_ANN_MODEL_H__
#define __DIRECT_ANN_MODEL_H__
#ifdef HAVE_CONFIG_H
#include "surfpack_config.h"
#endif

#include "surfpack_system_headers.h"
#include "SurfpackModel.h"
#include "SurfpackMatrix.h"


class DirectANNBasisSet
{
public:
  MtxDbl weights;
  DirectANNBasisSet(const MtxDbl& weights_in);
  double eval(unsigned index, const VecDbl& x) const;
  double deriv(unsigned index, const VecDbl& x, const VecUns& vars) const;
  double nodeSum(unsigned index, const VecDbl& x) const;
  std::string asString() const;
};

class DirectANNModel : public SurfpackModel
{
public:
  DirectANNModel(const DirectANNBasisSet& bs_in, const VecDbl& coeffs_in);
  virtual VecDbl gradient(const VecDbl& x) const;
  virtual std::string asString() const;
protected:
  virtual double evaluate(const VecDbl& x) const;
  DirectANNBasisSet bs;
  VecDbl coeffs;
friend class DirectANNModelTest;
};

///////////////////////////////////////////////////////////
///   Direct ANN Model Factory	
///////////////////////////////////////////////////////////

class DirectANNModelFactory : public SurfpackModelFactory 
{

public:
  DirectANNModelFactory();
  DirectANNModelFactory(const ParamMap& args);
  virtual SurfpackModel* Create(const SurfData& sd);
  virtual SurfpackModel* Create(const std::string& model_string);
  virtual void config();
protected:
  unsigned nNodes;
};
#endif
