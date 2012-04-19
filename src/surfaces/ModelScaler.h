/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#ifndef __MODEL_SCALER_H__
#define __MODEL_SCALER_H__

#include "surfpack_system_headers.h"

class SurfData;
class ModelScaler {
public:
  virtual const VecDbl& scale(const VecDbl& unscaled_x) const = 0;
  virtual double descale(double scaled_response) const = 0;
  virtual double scaleResponse(double unscaled_response) const = 0;
  virtual std::string asString() { return ""; }
  ModelScaler() {}
  virtual ~ModelScaler() {}
  virtual ModelScaler* clone() const = 0;
};

class NonScaler : public ModelScaler {
public:
  virtual const VecDbl& scale(const VecDbl& unscaled_x) const;
  virtual double descale(double scaled_response) const ;
  virtual double scaleResponse(double unscaled_response) const ;
  virtual std::string asString();
  static ModelScaler* Create(const SurfData& data);
  virtual ModelScaler* clone() const;
};

class NormalizingScaler : public ModelScaler {
public:
  struct Scaler {
    double offset;
    double scaleFactor;
    Scaler(double o, double s) : offset(o), scaleFactor(s) {}
    Scaler() {}
  };
  virtual const VecDbl& scale(const VecDbl& unscaled_x) const;
  virtual double descale(double scaled_response) const;
  virtual double scaleResponse(double unscaled_response) const ;
  virtual std::string asString();
  NormalizingScaler(const std::vector<Scaler>& s, const Scaler& d) 
    : scalers(s), descaler(d), result(s.size()) {}
  ~NormalizingScaler() {}
  // constructor to normalize each var/resp to [ 0, 1 ]
  static ModelScaler* Create(const SurfData& data);
  // constructor to normalize each var/resp to [ -norm_factor, norm_factor ]
  static ModelScaler* Create(const SurfData& data, double norm_factor);
  virtual ModelScaler* clone() const;
protected:
  std::vector<Scaler> scalers;
  Scaler descaler;
  mutable VecDbl result;
  friend class ModelScalerTest;
};

class ScaledSurfData {
public:
  ScaledSurfData(const ModelScaler& ms_in, const SurfData& sd_in);
  std::vector< double > getResponses() const;
  double getResponse(unsigned index) const;
  unsigned size() const;
  unsigned xSize() const;
  double operator()(unsigned pt, unsigned dim) const;
  const std::vector<double>& operator()(unsigned pt) const;
  static VecVecDbl asVecVecDbl(const ScaledSurfData& data);
protected:
  const ModelScaler& ms;
  const SurfData& sd;
};

#endif
