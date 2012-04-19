/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#ifndef __MARS_MODEL_H__
#define __MARS_MODEL_H__

#include "surfpack_system_headers.h"
#include "SurfpackModel.h"

typedef float real;


class MarsModel : public SurfpackModel
{
public:
  MarsModel(const unsigned dims, real* fm_in, int fmsize, int* im_in, 
    int imsize, int interp);
  virtual VecDbl gradient(const VecDbl& x) const;
  virtual std::string asString() const;
protected:
  virtual double evaluate(const VecDbl& x) const;
  VecDbl coeffs;

  std::vector<real> fm;
  std::vector<int> im;
  int interpolation;
friend class MarsModelTest;
};

///////////////////////////////////////////////////////////
///	Linear Regression Model Factory	
///////////////////////////////////////////////////////////

class MarsModelFactory : public SurfpackModelFactory 
{

public:
  MarsModelFactory();
  MarsModelFactory(const ParamMap& args);

protected:

  /// Model-specific portion of creation process
  virtual SurfpackModel* Create(const SurfData& sd);
  /// Model-specific portion of creation process
  virtual SurfpackModel* Create(const std::string& model_string);

  /// set member data prior to build; appeals to SurfpackModel::config()
  virtual void config();

  real* xMatrix;
  real* fm;
  int* im;
  int n;
  int np;
  int max_bases; 
  int max_interactions;
  int interpolation;
};
#endif
