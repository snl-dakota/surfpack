/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#ifndef __MOVING_LEAST_SQUARES_MODEL_H__
#define __MOVING_LEAST_SQUARES_MODEL_H__
#ifdef HAVE_CONFIG_H
#include "surfpack_config.h"
#endif

#include "surfpack_system_headers.h"
#include "SurfpackModel.h"
#include "SurfData.h"
#include "LinearRegressionModel.h"

class MovingLeastSquaresModel : public SurfpackModel
{
public:
  MovingLeastSquaresModel(const SurfData& sd_in, const LRMBasisSet& bs_in,
    unsigned continuity_in = 1);
  virtual VecDbl gradient(const VecDbl& x) const;
  virtual std::string asString() const;
protected:
  virtual double evaluate(const VecDbl& x) const;
  SurfData sd;
  LRMBasisSet bs;
  mutable VecDbl coeffs;
  unsigned continuity;
  
friend class MovingLeastSquaresModelTest;
};

///////////////////////////////////////////////////////////
///	Moving Least Squares Model Factory
///////////////////////////////////////////////////////////

class MovingLeastSquaresModelFactory : public SurfpackModelFactory 
{

public:
  MovingLeastSquaresModelFactory();
  MovingLeastSquaresModelFactory(const ParamMap& args);
  virtual SurfpackModel* Create(const SurfData& sd);
  virtual SurfpackModel* Create(const std::string& model_string);
  virtual void config();
protected:
  unsigned weight;
  unsigned order;
};
#endif
