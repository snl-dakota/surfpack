#include "config.h"

#ifndef SURF_SCALER_H
#define SURF_SCALER_H

#include <vector>

class SurfData;

class SurfScaler
{
public:
  /// To scale a value x, (compute x - offset) / divisor
  struct ScalingParameterPair {
    double offset;
    double divisor;
  };
  
  /// Initialize data to 0.;
  SurfScaler();
 
  /// Initialize data to data_, with scaledPoint set to appropriate size.
  SurfScaler(SurfData* data_);

  /// Initialize data to data_, with scaledPoint set to appropriate size.
  void setData(SurfData* data_);

  /// Iterate over a data set, computing the scaling parameters for each 
  /// dimension.
  void calculateParams();

  /// Return a scaled version of the parameter SurfPoint
  const std::vector<double>& operator()(const std::vector<double>& x) const;

private:
  std::vector<ScalingParameterPair> parameters;
  SurfData* data;
  mutable std::vector<double> scaledPoint;
   


};

#endif
