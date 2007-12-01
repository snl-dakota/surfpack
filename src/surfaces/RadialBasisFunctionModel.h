/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#ifndef __RADIAL_BASIS_FUNCTION_MODEL_H__
#define __RADIAL_BASIS_FUNCTION_MODEL_H__
#ifdef HAVE_CONFIG_H
#include "surfpack_config.h"
#endif

#include "surfpack_system_headers.h"
#include "SurfpackModel.h"
#include "SurfData.h"
#include "LinearRegressionModel.h"

class AxesBounds;
SurfPoint computeCentroid(const SurfData& sd);
void updateCentroid(VecDbl& centroid, const VecDbl& newpt, unsigned weight);
SurfData cvts(const AxesBounds& ab);
SurfData radii(const SurfData& generators);
VecUns probInclusion(unsigned vec_size, double prob);
VecDbl fullCoeff(unsigned vec_size, const VecDbl& coeffs, VecUns& incl);

class RadialBasisFunction
{
public:
  RadialBasisFunction(const VecDbl& center_in, const VecDbl& radius_in);
  RadialBasisFunction(const std::string& center_in, const std::string& radius_in);
  double operator()(const VecDbl& x) const;
  double deriv(const VecDbl& x, const VecUns& vars) const;
  std::string asString() const;
//protected:
  VecDbl center;
  VecDbl radius;
};
typedef std::vector<RadialBasisFunction> VecRbf;
VecRbf makeRbfs(const SurfData& generators, const SurfData& radii);
void augment(VecRbf& rbfs);


class RadialBasisFunctionModel : public SurfpackModel
{
public:
  RadialBasisFunctionModel(const VecRbf& rbfs_in, const VecDbl& coeffs_in);
  virtual double evaluate(const VecDbl& x) const;
  virtual VecDbl gradient(const VecDbl& x) const;
  virtual std::string asString() const;
protected:
  VecRbf rbfs;
  VecDbl coeffs;
friend class RadialBasisFunctionModelTest;
};

///////////////////////////////////////////////////////////
///	Radial Basis Function Model Factory
///////////////////////////////////////////////////////////

class RadialBasisFunctionModelFactory : public SurfpackModelFactory 
{

public:
  RadialBasisFunctionModelFactory();
  RadialBasisFunctionModelFactory(const ParamMap& args);
  virtual SurfpackModel* Create(const SurfData& sd);
  virtual SurfpackModel* Create(const std::string& model_string);
  virtual void config();
protected:
  unsigned nCenters;
  unsigned minPartition;
};

#endif
