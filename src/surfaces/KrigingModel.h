/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#ifndef __KRIGING_MODEL_H__
#define __KRIGING_MODEL_H__

#include "surfpack_system_headers.h"
#include "SurfpackModel.h"

#include "nkm/NKM_KrigingModel.hpp"
#include "nkm/NKM_GradKrigingModel.hpp"

class ScaledSurfData;

/// A thin wrapper around a NewKrigingModel
class KrigingModel : public SurfpackModel
{
public:
  KrigingModel(const SurfData& sd, const ParamMap& args);
  ~KrigingModel();
  virtual double variance(const VecDbl& x) const;
  virtual VecDbl gradient(const VecDbl& x) const;
  virtual MtxDbl hessian(const VecDbl& x) const;
  virtual std::string asString() const;
protected:
  MtxDbl getMatrix(const ScaledSurfData& ssd, const VecDbl& correlations);
  virtual double evaluate(const VecDbl& x) const;
friend class KrigingModelTest;

private:
  // helper
  void surfdata_to_nkm_surfdata(const SurfData& sd, nkm::SurfData& nkm_sd);

  // use class data to keep in scope for wrapped model
  //nkm::KrigingModel* nkmKrigingModel;
  nkm::SurfPackModel* nkmKrigingModel;
  nkm::SurfData nkmSurfData;
};


///////////////////////////////////////////////////////////
///   Kriging Model Factory	
///////////////////////////////////////////////////////////

class KrigingModelFactory : public SurfpackModelFactory 
{

public:
  KrigingModelFactory();
  KrigingModelFactory(const ParamMap& args);
  /// Override since Kriging does allow constraints
  virtual bool supports_constraints();

protected:

  /// Model-specific portion of creation process
  virtual SurfpackModel* Create(const SurfData& sd);
  /// Model-specific portion of creation process
  virtual SurfpackModel* Create(const std::string& model_string);

  /// set member data prior to build; appeals to SurfpackModel::config()
  virtual void config();

  /// For Kriging, sufficient data is assessed by the NKM submodel
  virtual void sufficient_data(const SurfData& sd);
};
#endif
