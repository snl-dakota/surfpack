/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#ifndef __DIRECT_ANN_MODEL_H__
#define __DIRECT_ANN_MODEL_H__

#include "surfpack_system_headers.h"
#include "SurfpackModel.h"

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
  MtxDbl randomMatrix(unsigned nrows, unsigned ncols);

protected:

  /// Model-specific portion of creation process
  virtual SurfpackModel* Create(const SurfData& sd);
  /// Model-specific portion of creation process
  virtual SurfpackModel* Create(const std::string& model_string);

  /// set member data prior to build; appeals to SurfpackModel::config()
  virtual void config();

  unsigned nodes;
  double range;
  unsigned samples;
};
#endif
