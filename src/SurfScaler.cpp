#include "config.h"

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <vector>

#include "surfpack.h"
#include "SurfPoint.h"
#include "SurfData.h"
#include "SurfScaler.h"

using namespace std;
SurfScaler::SurfScaler()
  : scaledPoint(vector<double>(1)),
  parameters(0)
{

}

SurfScaler::SurfScaler(const SurfScaler& other)
  : scaledPoint(other.scaledPoint), 
  parameters(other.parameters)
{

}

void SurfScaler::computeScalingParameters(const SurfData& sd)
{
  if (sd.size() > 0) {
    scaledPoint.resize(sd.xSize());
    parameters.resize(sd.xSize());
    for (unsigned i = 0; i < parameters.size(); i++ ) {
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
      parameters[i].offset = min;
      // the divisor is the range of the values along this dimension
      parameters[i].divisor = max - min;
    }
  }
}

/// Return a scaled version of the parameter SurfPoint
const SurfPoint& SurfScaler::scale(const vector<double>& x) const
{
  if ((x.size() == parameters.size() && 
       x.size() == scaledPoint.xSize())) {
    for (unsigned i = 0; i < x.size(); i++) {
      scaledPoint.setX(i,
        (x[i] - parameters[i].offset) / parameters[i].divisor);
    }
  }
  return scaledPoint;
}
