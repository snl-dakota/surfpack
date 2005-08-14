#include "surfpack_config.h"
#include "surfpack.h"
#include "SurfData.h"
#include "SurfScaler.h"

using namespace std;

DimensionScaler::DimensionScaler()
{

}

DimensionScaler::~DimensionScaler()
{

}

NonScaler::NonScaler()
{

}

NonScaler::~NonScaler()
{

}

double NonScaler::scale(double value)
{
  return value;
}

double NonScaler::descale(double value)
{
  return value;
}

DimensionScaler* NonScaler::clone()
{
  return new NonScaler();
}

std::string NonScaler::asString()
{
  return string("Non-scaler");
}
NormalizingScaler::NormalizingScaler(const NormalizingScaler& other)
  : offset(other.offset), divisor(other.divisor)
{

}

NormalizingScaler::NormalizingScaler(const SurfData& surf_data, 
  unsigned dim, bool response) : offset(0.0), divisor(1.0)
{
  if (surf_data.size() > 0) {
    double vmin;
    double vmax;
    if (!response) {
      assert(dim < surf_data.xSize());
      // find the vmax and vmin values along ith dimension over all j data 
      // points.  sd.scaler should be null until after this method
      // returns.
      vmin = surf_data[0][dim];
      vmax = surf_data[0][dim];
      for (unsigned j = 1; j < surf_data.size(); j++ ) {
        if (surf_data[j][dim] > vmax) {
          vmax = surf_data[j][dim];
        } else if (surf_data[j][dim] < vmin) {
          vmin = surf_data[j][dim];
        }
      }
    } else {
      assert(dim < surf_data.fSize());
      // find the vmax and vmin values along ith dimension over all j data 
      // points.  sd.scaler should be null until after this method
      // returns.
      vmin = surf_data[0].F(dim);
      vmax = surf_data[0].F(dim);
      for (unsigned j = 1; j < surf_data.size(); j++ ) {
        if (surf_data[j].F(dim) > vmax) {
          vmax = surf_data[j].F(dim);
        } else if (surf_data[j].F(dim) < vmin) {
          vmin = surf_data[j].F(dim);
        }
      }
    }
    // the offset for this dimension is the minimum value
    offset = vmin;
    // the divisor is the range of the values along this dimension
    divisor = vmax - vmin;
    //cout << "For design var " << i << " offset: " << designVarParams[i].offset
    //     << " divisor: " << designVarParams[i].divisor << endl;
  }
}

NormalizingScaler::~NormalizingScaler()
{

}

double NormalizingScaler::scale(double value)
{
  return (value - offset) / divisor;
}

double NormalizingScaler::descale(double value)
{
  return (value * divisor) + offset;
}

DimensionScaler* NormalizingScaler::clone()
{
  return new NormalizingScaler(*this);
}

std::string NormalizingScaler::asString()
{
  ostringstream os;
  os << "offset: " << offset << " divisor: " << divisor;
  return os.str();
}

SurfScaler::SurfScaler()
{

}

SurfScaler::SurfScaler(const SurfScaler& other)
{
  scalers.resize(other.scalers.size());
  responseScalers.resize(other.responseScalers.size());
  for(unsigned i = 0; i < scalers.size(); i++) {
    scalers[i] = other.scalers[i]->clone();
  }
  for(unsigned i = 0; i < responseScalers.size(); i++) {
    responseScalers[i] = other.responseScalers[i]->clone();
  }
}

/// Delete members of scalers and responseScalers
SurfScaler::~SurfScaler()
{
  for(unsigned i = 0; i < scalers.size(); i++) {
    delete scalers[i];
    scalers[i] = 0;
  }
  for(unsigned i = 0; i < responseScalers.size(); i++) {
    delete responseScalers[i];
    responseScalers[i] = 0;
  }
}
 
void SurfScaler::computeScalingParameters(const SurfData& sd)
{
  scalers = vector< DimensionScaler* >(sd.xSize(),0);
  responseScalers = vector< DimensionScaler* >(sd.fSize(),0);
  if (sd.size() > 0) {
    for(unsigned i = 0; i < scalers.size(); i++) {
      scalers[i] = new NormalizingScaler(sd,i);
    }
    for(unsigned i = 0; i < responseScalers.size(); i++) {
      responseScalers[i] = new NormalizingScaler(sd,i,true);
    }
  }
}

/// Return a scaled version of the parameter SurfPoint
double SurfScaler::scale(unsigned index, double value) const
{
  assert(index < scalers.size());
  return scalers[index]->scale(value);
}

/// Return the scaled value of a response variable
double SurfScaler::scaleResponse(unsigned index, double value)
{
  assert(index < responseScalers.size());
  return responseScalers[index]->scale(value);
}

double SurfScaler::descaleResponse(unsigned index, double value)
{
  assert(index < responseScalers.size());
  return responseScalers[index]->descale(value);
}

std::string SurfScaler::asString()
{
  ostringstream os;
  os << "predictor variables" << endl;
  for (unsigned i = 0; i < scalers.size(); i++) {
    os << i << " " << scalers[i]->asString() << endl;
  }
  os << "response variables" << endl;
  for (unsigned i = 0; i < responseScalers.size(); i++) {
    os << i << " " << responseScalers[i]->asString() << endl;
  }
  return os.str();
}
