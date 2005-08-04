#include "surfpack_config.h"
#include "surfpack.h"
#include "SurfData.h"
#include "SurfScaler.h"

using namespace std;
SurfScaler::SurfScaler()
  : scaledPoint(vector<double>(1)),
  designVarParams(0), responseVarParams(0)
{

}

SurfScaler::SurfScaler(const SurfScaler& other)
  : scaledPoint(other.scaledPoint), 
  designVarParams(other.designVarParams),
  responseVarParams(other.responseVarParams)
{

}

void SurfScaler::computeScalingParameters(const SurfData& sd)
{
  if (sd.size() > 0) {
    scaledPoint.resize(sd.xSize());
    designVarParams.resize(sd.xSize());
    for (unsigned i = 0; i < designVarParams.size(); i++ ) {
      // find the max and min values along ith dimension over all j data 
      // points.  sd.scaler should be null until after this method
      // returns.
      double min = sd[0].X()[i];
      double max = sd[0].X()[i];
      for (unsigned j = 1; j < sd.size(); j++ ) {
        if (sd[j].X()[i] > max) {
          max = sd[j].X()[i];
        } else if (sd[j].X()[i] < min) {
          min = sd[j].X()[i];
        }
      }
      // the offset for this dimension is the minimum value
      designVarParams[i].offset = min;
      // the divisor is the range of the values along this dimension
      designVarParams[i].divisor = max - min;
      //cout << "For design var " << i << " offset: " << designVarParams[i].offset
      //     << " divisor: " << designVarParams[i].divisor << endl;
    }

    responseVarParams.resize(sd.fSize());
    // Repeat the same process for the response variables
    for (unsigned i = 0; i < responseVarParams.size(); i++ ) {
      // find the max and min values along ith dimension over all j data 
      // points.  sd.scaler should be null until after this method
      // returns.
      double min = sd[0].F(i);
      double max = sd[0].F(i);
      for (unsigned j = 1; j < sd.size(); j++ ) {
        if (sd[j].F(i) > max) {
          max = sd[j].F(i);
        } else if (sd[j].F(i) < min) {
          min = sd[j].F(i);
        }
      }
      // the offset for this dimension is the minimum value
      responseVarParams[i].offset = min;
      // the divisor is the range of the values along this dimension
      responseVarParams[i].divisor = max - min;
      //cout << "For response var " << i << " offset: " << responseVarParams[i].offset
      //     << " divisor: " << responseVarParams[i].divisor << endl;
    }
  }
}

/// Return a scaled version of the parameter SurfPoint
const SurfPoint& SurfScaler::scale(const vector<double>& x) const
{
  if ((x.size() == designVarParams.size() && 
       x.size() == scaledPoint.xSize())) {
    for (unsigned i = 0; i < x.size(); i++) {
      scaledPoint.setX(i,
        (x[i] - designVarParams[i].offset) / designVarParams[i].divisor);
    }
  }
  return scaledPoint;
}

/// Return the scaled value of a response variable
double SurfScaler::scaleResponse(double value, unsigned index)
{
  if (index < 0 || index >= responseVarParams.size()) {
    throw string(
      "Range error in scaling: index out of range in SurfScaler::scaleResponse"
    );
  }
  return (value - responseVarParams[index].offset) / 
    responseVarParams[index].divisor;
}

double SurfScaler::descaleResponse(double value, unsigned index)
{
  if (index < 0 || index >= responseVarParams.size()) {
    throw string(
      "Range error in descaling: index out of range in descaleResponse"
    );
  }
  return (value * responseVarParams[index].divisor) + 
    responseVarParams[index].offset;
}
