/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#ifndef __KRIGING_MODEL_H__
#define __KRIGING_MODEL_H__
#ifdef HAVE_CONFIG_H
#include "surfpack_config.h"
#endif

#include "surfpack_system_headers.h"
#include "SurfpackModel.h"
#include "SurfpackMatrix.h"


class KrigingBasisSet
{
public:
  VecVecDbl centers;
  VecDbl correlations;
  KrigingBasisSet(const VecVecDbl& centers_in, const VecDbl& correlations_in);
  double eval(unsigned index, const VecDbl& x) const;
  double deriv(unsigned index, const VecDbl& x, const VecUns& vars) const;
  std::string asString() const;
};

class KrigingModel : public SurfpackModel
{
public:
  KrigingModel(const KrigingBasisSet& bs_in, const VecDbl& rhs_in);
  KrigingModel(const SurfData& sd, const VecDbl& correlations);
  virtual VecDbl gradient(const VecDbl& x) const;
  virtual std::string asString() const;

  double likelihood;
  static MtxDbl corrMtx(const VecDbl& corr_vec, const SurfData& data);
protected:
  MtxDbl getMatrix(const ScaledSurfData& ssd, const VecDbl& correlations);
  virtual double evaluate(const VecDbl& x) const;
  KrigingBasisSet bs;
  VecDbl rhs;
  double betaHat;
friend class KrigingModelTest;
};

#include "Conmin.h"
class ConminKriging : public Conmin
{
public:
  ConminKriging(const SurfData& data_in);
  virtual void optimize(VecDbl& x, double& final_val, unsigned max_iter);
  virtual double objective(const VecDbl& x);
  virtual VecDbl gradient(const VecDbl& x);
protected:
  const SurfData& data;
  VecDbl rhs;
};

///////////////////////////////////////////////////////////
///   Kriging Model Factory	
///////////////////////////////////////////////////////////

class KrigingModelFactory : public SurfpackModelFactory 
{

public:
  KrigingModelFactory();
  KrigingModelFactory(const ParamMap& args);
  virtual SurfpackModel* Create(const SurfData& sd);
  virtual SurfpackModel* Create(const std::string& model_string);
  virtual void config();
protected:
  VecDbl sampleCorrelations(const SurfData& sd);
  VecDbl conminCorrelations(const SurfData& sd);
  VecDbl correlations;
  int max_iter;
  VecDbl conmin_seed;
};
#endif
