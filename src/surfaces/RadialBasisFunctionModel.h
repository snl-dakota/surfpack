/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#ifndef __RADIAL_BASIS_FUNCTION_MODEL_H__
#define __RADIAL_BASIS_FUNCTION_MODEL_H__

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
  /// default constructor used when reading from archive file 
  RadialBasisFunction() { /* empty ctor */}
  RadialBasisFunction(const VecDbl& center_in, const VecDbl& radius_in);
  RadialBasisFunction(const std::string& center_in, const std::string& radius_in);
  double operator()(const VecDbl& x) const;
  double deriv(const VecDbl& x, const VecUns& vars) const;
  std::string asString() const;
//protected:
  VecDbl center;
  VecDbl radius;

private:

  // allow serializers access to private data
  friend class boost::serialization::access;
  /// serializer for derived class Model data
  template<class Archive> 
  void serialize(Archive & archive, const unsigned int version);

};
typedef std::vector<RadialBasisFunction> VecRbf;
VecRbf makeRbfs(const SurfData& generators, const SurfData& radii);
void augment(VecRbf& rbfs);


class RadialBasisFunctionModel : public SurfpackModel
{
public:
  RadialBasisFunctionModel() { /* empty ctor */ }
  RadialBasisFunctionModel(const VecRbf& rbfs_in, const VecDbl& coeffs_in);
  virtual double evaluate(const VecDbl& x) const;
  virtual VecDbl gradient(const VecDbl& x) const;
  virtual std::string asString() const;
protected:
  VecRbf rbfs;
  VecDbl coeffs;
friend class RadialBasisFunctionModelTest;

private:

  // allow serializers access to private data
  friend class boost::serialization::access;
  /// serializer for derived class Model data
  template<class Archive> 
  void serialize(Archive & archive, const unsigned int version);

};

BOOST_CLASS_EXPORT_KEY(RadialBasisFunctionModel)

///////////////////////////////////////////////////////////
///	Radial Basis Function Model Factory
///////////////////////////////////////////////////////////

class RadialBasisFunctionModelFactory : public SurfpackModelFactory 
{

public:
  RadialBasisFunctionModelFactory();
  RadialBasisFunctionModelFactory(const ParamMap& args);
protected:

  /// Model-specific portion of creation process
  virtual SurfpackModel* Create(const SurfData& sd);
  /// Model-specific portion of creation process
  virtual SurfpackModel* Create(const std::string& model_string);

  /// set member data prior to build; appeals to SurfpackModel::config()
  virtual void config();

  unsigned nCenters;
  unsigned cvtPts;
  unsigned maxSubsets;
  unsigned minPartition;
};

#endif
