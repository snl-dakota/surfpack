/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
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
  MtxDbl getMatrix(const ScaledSurfData& ssd, const VecDbl& correlations);
  virtual VecDbl gradient(const VecDbl& x) const;
  virtual std::string asString() const;
  static KrigingModel Create(const SurfData& sd);
  static KrigingModel Create(const SurfData& sd, const VecDbl& correlations);
  static MtxDbl corrMtx(const VecDbl& corr_vec, const SurfData& data);
protected:
  virtual double evaluate(const VecDbl& x) const;
  KrigingBasisSet bs;
  VecDbl rhs;
  double betaHat;
  double likelihood;
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
#endif
