#ifndef SURF_SCALER_H
#define SURF_SCALER_H

#include "surfpack_config.h"
#include "surfpack_system_headers.h"

class SurfData;
#include "SurfPoint.h"

class DimensionScaler
{
public:
  DimensionScaler();
  virtual ~DimensionScaler();
  virtual double scale(double value) = 0;
  virtual double descale(double value) = 0;
  virtual DimensionScaler* clone() = 0;
  virtual std::string asString() = 0;
#ifdef __TESTING_MODE__
  friend class SurfScalerTest;
  friend class SurfDataTest;
  friend class SurfaceTest;
#endif
};

class NonScaler : public DimensionScaler
{
public:
  NonScaler();
  virtual ~NonScaler();
  virtual double scale(double value);
  virtual double descale(double value);
  virtual DimensionScaler* clone();
  virtual std::string asString();
};

class NormalizingScaler : public DimensionScaler
{
public:
  NormalizingScaler(const SurfData& surf_data, unsigned dim, 
    bool response = false);
  NormalizingScaler(const NormalizingScaler& other);
  virtual ~NormalizingScaler();
  virtual double scale(double value);
  virtual double descale(double value);
  virtual DimensionScaler* clone();
  virtual std::string asString();
#ifdef __TESTING_MODE__
  friend class SurfScalerTest;
  friend class SurfDataTest;
  friend class SurfaceTest;
#endif
protected:
  /// To scale a value x, (compute x - offset) / divisor
  double offset;
  double divisor;
};

class SurfScaler
{
public:
  /// Initializes scaledPoint 
  SurfScaler();

  /// Makes deep copy
  SurfScaler(const SurfScaler& other);

  /// Delete members of scalers and responseScalers
  ~SurfScaler();
 
  /// Iterate over a data set, computing the scaling parameters for each 
  /// dimension.
  void computeScalingParameters(const SurfData& data_);

  /// Return a scaled version of the parameter point 
  double scale(unsigned index, double value) const;

  ///// Return the scaled value of a response variable
  double scaleResponse(unsigned index, double value);

  ///// Invert the scaling on the response variable
  double descaleResponse(unsigned index, double value);

  /// Express SurfScaler as a string, with one DimensionScaler per line
  std::string asString();

  

protected:
  std::vector< DimensionScaler* > scalers;
  std::vector< DimensionScaler* > responseScalers;
   
#ifdef __TESTING_MODE__
  friend class SurfScalerTest;
  friend class SurfDataTest;
  friend class SurfaceTest;
#endif
};

#endif
