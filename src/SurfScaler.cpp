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
  : data(0), scaledPoint()
{

}

SurfScaler::SurfScaler(SurfData* data_)
  : data(0), scaledPoint()
{
  setData(data_);
}

void SurfScaler::setData(SurfData* data_)
{
  this->data = data_;
  // initialize the scaled point to have the right number of dimensions
  scaledPoint.resize(data->size());
}

void SurfScaler::calculateParams()
{
  if (data) {
    SurfData& sd = *data;
    if (sd.size() > 0) {
      parameters.resize(sd.xSize());
      for (unsigned i = 0; i < parameters.size(); i++ ) {
        // find the max and min values along ith dimension over all j data 
        // points.
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
}

/// Return a scaled version of the parameter SurfPoint
const vector<double>& SurfScaler::operator()
  (const vector<double>& x) const
{
  if ((x.size() == parameters.size() && 
       x.size() == scaledPoint.size())) {
    for (unsigned i = 0; i < x.size(); i++) {
      scaledPoint[i] = (x[i] - parameters[i].offset) / parameters[i].divisor;
    }
  }
  return scaledPoint;
}
