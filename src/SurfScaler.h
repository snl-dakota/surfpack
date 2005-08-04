#ifndef SURF_SCALER_H
#define SURF_SCALER_H

#include "surfpack_config.h"
#include "surfpack_system_headers.h"

class SurfData;
#include "SurfPoint.h"

class SurfScaler
{
public:
  /// To scale a value x, (compute x - offset) / divisor
  struct ScalingParameterPair {
    double offset;
    double divisor;
  };
  
  /// Initializes scaledPoint 
  SurfScaler();

  /// Makes deep copy
  SurfScaler(const SurfScaler& other);
 
  /// Iterate over a data set, computing the scaling parameters for each 
  /// dimension.
  void computeScalingParameters(const SurfData& data_);

  /// Return a scaled version of the parameter point 
  const SurfPoint& scale(const std::vector<double>& x) const;

  /// Return the scaled value of a response variable
  double scaleResponse(double value, unsigned index);

  /// Invert the scaling on the response variable
  double descaleResponse(double value, unsigned index);

protected:
  std::vector<ScalingParameterPair> designVarParams;
  std::vector<ScalingParameterPair> responseVarParams;
  mutable SurfPoint scaledPoint;
   
#ifdef __TESTING_MODE__
  friend class SurfScalerTest;
  friend class SurfDataTest;
  friend class SurfaceTest;
#endif


};

#endif
