#include "config.h"

#ifndef SURF_SCALER_H
#define SURF_SCALER_H

#include <vector>

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

protected:
  std::vector<ScalingParameterPair> parameters;
  mutable SurfPoint scaledPoint;
   
#ifdef __TESTING_MODE__
  friend class SurfScalerTest;
  friend class SurfDataTest;
  friend class SurfaceTest;
#endif


};

#endif
