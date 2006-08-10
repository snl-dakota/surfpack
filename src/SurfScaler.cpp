/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

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

DimensionScaler* NonScaler::clone() const
{
  return new NonScaler();
}

string NonScaler::asString()
{
  return string("Non-scaler");
}

LogScaler::LogScaler()
{

}

LogScaler::~LogScaler()
{

}

double LogScaler::scale(double value)
{
  return log(value);
}

double LogScaler::descale(double value)
{
  return exp(value);
}

DimensionScaler* LogScaler::clone() const
{
  return new LogScaler();
}

string LogScaler::asString()
{
  return string("Log-scaler");
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

DimensionScaler* NormalizingScaler::clone() const
{
  return new NormalizingScaler(*this);
}

string NormalizingScaler::asString()
{
  ostringstream os;
  os << "offset: " << offset << " divisor: " << divisor;
  return os.str();
}

// ____________________________________________________________________________
// Creation, Destruction, Initialization 
// ____________________________________________________________________________

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
  cleanup();
}

/// Delete members of scalers and responseScalers
void SurfScaler::cleanup()
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
 
// ____________________________________________________________________________
// Queries 
// ____________________________________________________________________________

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

string SurfScaler::asString()
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

/// Number of predictor variables that this object scales
unsigned SurfScaler::xSize()
{
  return scalers.size();
}

/// Number of responses that this object scales
unsigned SurfScaler::fSize()
{
  return responseScalers.size();
}

// ____________________________________________________________________________
// Commands 
// ____________________________________________________________________________

void SurfScaler::normalizeAll(const SurfData& sd)
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

/// Sets the scaler for a response or predictor variable.  
void SurfScaler::setDimensionScaler(unsigned index, const DimensionScaler& ds, 
  bool response)
{
  if (response) {
    assert(index < responseScalers.size());
    delete responseScalers[index];
    responseScalers[index] = ds.clone();
  } else {
    assert(index < scalers.size());
    delete scalers[index];
    scalers[index] = ds.clone();
  }
}

/// Use the same DimensionScaler for all variables
void SurfScaler::setAll(const DimensionScaler& ds)
{
    for(unsigned i = 0; i < scalers.size(); i++) {
      delete scalers[i];
      scalers[i] = ds.clone();
    }
    for(unsigned i = 0; i < responseScalers.size(); i++) {
      delete responseScalers[i];
      responseScalers[i] = ds.clone();
    }
}

/// Let Surfpack guess how to best scale each variable
void SurfScaler::scaleAuto(const SurfData& sd)
{
  normalizeAll(sd);
}

/// Use SurfpackParserArgs to determine scaling for each parameter
void SurfScaler::configList(const SurfData& sd, ArgList& args)
{
  for (unsigned i = 0; i < args.size(); i++) {
    config(sd,args[i]);
  }
}

/// Parse out info from arg to configure SurfScaler
void SurfScaler::config(const SurfData& sd, const Arg& arg)
{
  bool response;
  unsigned index;
  sizeToMatch(sd);
  if (arg.name == "norm_scale") {
    //\todo Check for arglist
    // If yes, parse out offset and divisor
    // parse out list of args 
    const Tuple& vars = arg.getRVal()->getTuple();
    for (unsigned i = 0; i < vars.size(); i++) {
      bool found = sd.varIndex(vars[i],index,response);
      assert(found);
      if (found) {
        setDimensionScaler(index, NormalizingScaler(sd,index,response), 
          response);
      }
    }
  }
  if (arg.name == "log_scale") {
    const Tuple& vars = arg.getRVal()->getTuple();
    for (unsigned i = 0; i < vars.size(); i++) {
      bool found = sd.varIndex(vars[i],index,response);
      assert(found);
      if (found) {
        setDimensionScaler(index, LogScaler(), 
          response);
      }
    }
  }
}

/// Size scalers and responseScalers to match data set; initialize
/// each dim to non-scaling
void SurfScaler::sizeToMatch(const SurfData& sd)
{
  NonScaler ns;
  while (scalers.size() != sd.xSize()) {
    scalers.push_back(ns.clone());
  }
  while (responseScalers.size() != sd.fSize()) {
    responseScalers.push_back(ns.clone());
  }
}
  
// ____________________________________________________________________________
// ____________________________________________________________________________
// Helper methods 
// ____________________________________________________________________________

