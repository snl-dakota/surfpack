#include "config.h"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "surfpack.h"
#include "SurfPoint.h"
#include "SurfData.h"
#include "Surface.h"
#include "SurfScaler.h"

using namespace std;

// ____________________________________________________________________________
// Creation, Destruction, Initialization 
// ____________________________________________________________________________

/// Data to be used to create surface is specified 
// sd must be initialized to NULL before setData(...) is called, because 
// setData(...) will call delete on the old value of sd, which will be
// garbage if it's not initialized to NULL.
Surface::Surface(SurfData* sd_) 
  : sd(0), scaler(0)
{ 
  setData(sd_);
  init();
}

/// Deep copies are made of most data members.  For the sd member, the new
/// object simply asks to be added to sd's list of listeners. 
// sd must be initialized to NULL before setData(...) is called, because 
// setData(...) will call delete on the old value of sd, which will be
// garbage if it's not initialized to NULL.
Surface::Surface(const Surface& other) 
  : sd(0), scaler(other.scaler), xsize(other.xsize), 
  builtOK(other.builtOK), dataModified(other.dataModified), 
  responseIndex(other.responseIndex)
{
  cout << "Surface copy constructor called" << endl;
  setData(other.sd);
}

/// Notifies its SurfData object that it is no longer observing 
Surface::~Surface() 
{ 
  if(sd) {
     sd->removeListener(this);
  } 
  delete scaler;
  scaler = 0;
}

/// Common initialization for new objects 
void Surface::init()
{
  dataModified = builtOK = false;
  xsize = sd ? sd->xSize() : 0;
  responseIndex = 0;
}

// ____________________________________________________________________________
// Queries 
// ____________________________________________________________________________

/// Return dimensionality of the surface or zero if not built
unsigned Surface::xSize()
{
  return xsize;
}

/// Return true if the data that was used to create the surface is available.
/// Some error metrics require the original data.
bool Surface::hasOriginalData() const
{
  return (sd && builtOK && !dataModified); 
}

/// Return true if there is a sufficient number of correctly formatted data
/// points
bool Surface::acceptableData() const
{
  if (!sd) {
    throw SurfData::bad_surf_data("Data unacceptable: there is no data.");
  } else {
    unsigned pointsAvailable = sd->size();
    unsigned pointsRequired = minPointsRequired();
    if (pointsAvailable < pointsRequired) {
      ostringstream errormsg;
      errormsg << "ERROR: data unacceptable.  This surface requires "
	   << pointsRequired << ", but only " << pointsAvailable
	   << " were given." << endl;
      throw SurfData::bad_surf_data(errormsg.str());
    }
  }
  // If no exception has been thrown when execution gets here, the data is ok.
  return true;
}
  
/// Evaluate the approximation surface at point x and return the value.
/// The point x must have the same dimensionality as this surface's SurfData.
/// Make sure the surface is valid and then call evaluate
double Surface::getValue(const std::vector<double>& x)
{
  if (!builtOK || dataModified) {
    createModel();
  }
  if (!scaler) {
    return evaluate(x);
  } else {
    return scaler->descaleResponse(
      evaluate(scaler->scale(x).X()),responseIndex
    );
  }
}

/// Evaluate the approximation surface at point x and return the value.
/// The point x must have the same dimensionality as this surface's SurfData.
double Surface::getValue(const SurfPoint& sp)
{
  return getValue(sp.X());
}

/// Evaluate the approximation surface at each point in the parameter
/// SurfData object.  In the ErrorStruct list, store the expected value (as
/// returned by sd.getResponse()) and the estimated value.
void Surface::getValue(SurfData& surf_data, std::vector<surfpack::ErrorStruct>& pts)
{
  for (unsigned i = 0; i < surf_data.size(); i++) {
    surfpack::ErrorStruct es;
    es.estimated = getValue(surf_data[i].X());
    es.observed = surf_data.getResponse(i);
    pts.push_back(es);
  }
}

/// Evaluate the approximation surface at each point in the parameter
/// surfData object.  Append the evaluations as a new response variable in
/// the data set.
void Surface::getValue(SurfData& surfData)
{
  set<unsigned> emptySet;
  surfData.setExcludedPoints(emptySet);
  vector<double> newValues;
  for (unsigned i = 0; i < surfData.size(); i++) {
    newValues.push_back(getValue(surfData[i].X()));
  }
  surfData.addResponse(newValues);
}

/// Modify one of the configuration options for this surface type 
void Surface::config(const SurfpackParser::Arg& arg)
{
  // Other arguments are processed by config methods in the child class.
  // If the child class sees an argument it does not recognize, it should
  // call this method.
  if (arg.name == "response_index") {
    unsigned index = static_cast<unsigned>(arg.rval.integer);
    if (!sd) {
      throw string("Cannot set response index on NULL data");
    } else {
      sd->setDefaultIndex(index);
      this->responseIndex = index; 
    }
  } else if (arg.name == "scaling") {
    string scaleOption = arg.rval.literal;
    if (scaleOption == "uniform") {
      scaleUniform();
      cout << "Uniform scaling" << endl;
    } else if (scaleOption == "none") {
      noScale();
    } else {
      throw string("Unrecognized option for surface parameter 'scaling'");
    }  
  } else if (arg.name == "xsize") {
    setXSize(arg.rval.integer);
  }
}

/// Process a list of configuration options
void Surface::configList(const SurfpackParser::ArgList& arglist)
{
  for (unsigned i = 0; i < arglist.size(); i++) {
    config(arglist[i]);
  }
}

/// Sets the number of dimensions in the surface
void Surface::setXSize(unsigned xsize_in)
{
  this->xsize = xsize_in;
}

/// Return the value of some error metric
double Surface::goodnessOfFit(const string metricName, SurfData* surfData)
{
  SurfData& sdRef = checkData(surfData);
  if (metricName == "sse") {
    return sse(sdRef);
  } else if (metricName == "mse") {
    return mse(sdRef);
  } else if (metricName == "mrae") {
    return mrae(sdRef);
  } else if (metricName == "rsquared") {
    return rSquared(sdRef);
  } else if (metricName == "press") {
    return press(sdRef);
  } else {
    throw string("No error metric of that type in this class");
  }
}

/// For each point x in dataSet, construct a similar approximation Surface 
/// that includes all of the points in dataSet except x.  Then evaluate the
/// Surface at x.  The difference between x and the estimate of x given the
/// rest of the data is the residual for x.  The PRESS statistic is the 
/// square root of the mean of the squares of all the residuals.
double Surface::press(SurfData& dataSet)
{
  /// <= test is used because it must be possible to build the surface
  /// even when one point is removed from dataSet
  if (dataSet.size() <= minPointsRequired()) {
    throw SurfData::bad_surf_data("Not enough data to compute PRESS.");
  } else {
    // If some of the points in the data set are already being excluded,
    // copy all of the non-excluded data points into a new SurfData
    // object where none of the points are excluded.  This will make
    // it easier to do the "leave one out" process
    bool containsInactives = !dataSet.getExcludedPoints().empty();
    SurfData activeSet = (containsInactives) 
      ? dataSet.copyActive() : dataSet;
    double pressValue = 0.0;
    unsigned i = 0;

    // For each data point, make a new surface that has all of the
    // data points except the current one.  Evaluate the new
    // surface at the omitted point and compute the residual
    // between the actual value and the value predicted by the new 
    // model at That point.
    unsigned totalPoints = activeSet.size();
    i = 0;
    set<unsigned> pointToSkip;
    while (i < totalPoints) {
      pointToSkip.clear();
      activeSet.setExcludedPoints(pointToSkip);
      // Get the point from the full data set corresponding to the index that 
      // will be skipped.  (This point will not be available to test after it
      // has been marked for exclusion).
      const SurfPoint& currentPoint = activeSet[i];
      pointToSkip.insert(i);
      activeSet.setExcludedPoints(pointToSkip);

      Surface* allButOne = makeSimilarWithNewData(&activeSet);
      double fTilde = allButOne->getValue(currentPoint.X());
      // actual value of the function at this point
      double fx = currentPoint.F(this->responseIndex);
      // add the square of the residual (fTilde - fx)
      pressValue += (fTilde - fx) * (fTilde - fx);
      delete allButOne;
      allButOne = 0;
      i++;
    }
    pressValue = sqrt(pressValue/static_cast<double>(activeSet.size()+1));
    return pressValue;
  }
  return 0.0;
}

/// Statistically speaking, R^2 is extra sum of squares divided by the total
/// sum of squares.  It measures how much of the variation in the data is 
/// accounted for by the model (the approximating surface).
double Surface::rSquared(SurfData& dataSet)
{
  // Sum of the function evaluations for all the points; used to compute mean
  double sumObserved = 0.0;
  
  // Sum of the squared function evaluations for the points; used to compute
  // the total sum of squares
  double sumOfSquaresObserved = 0.0;

  // Sum of the squared differences between the true function value and the
  // Surface approximation's estimate of the function value over all of the
  // points
  double residualSumOfSquares = 0.0;

  // Sum of the squared differences between the true function value and mean 
  // function value over all of the points
  double totalSumOfSquares = 0.0;

  for (unsigned i = 0; i < dataSet.size(); i++) {
    double observedF = dataSet.getResponse(i);
    double estimatedF = getValue(dataSet[i].X());
    double residual = observedF - estimatedF;
    residualSumOfSquares += residual * residual;
    sumObserved += observedF;
    sumOfSquaresObserved += observedF * observedF;
  }
  // This is the same as sigma{i=1 to n}(x_i - xbar)^2
  totalSumOfSquares = sumOfSquaresObserved - 
    (sumObserved * sumObserved / dataSet.size());
  double rSquaredValue = 1.0 - residualSumOfSquares / totalSumOfSquares;
  // In a polynomial regression, this computation will always result in a value
  // between 0 and 1, because the residual sum of squares will always be less
  // than or equal to the total sum of squares (i.e., the worst your regression
  // can do is give you the mean everywhere).  Some of the other methods (Mars,
  // Kriging, RBF, etc.) can have much larger residual sum of squares.  We will
  // refrain from returning a value less than zero however.
  return (rSquaredValue < 0) ? 0 : rSquaredValue;
} 
      
/// The sum of squared errors.  The response variable at dataSet's 
/// defaultIndex is interpreted as the true function value.
double Surface::sse(SurfData& dataSet)
{
  vector<surfpack::ErrorStruct> results;
  getValue(dataSet,results);
  double sse = 0.0;
  for (unsigned i = 0; i < results.size(); i++) {
    double residual = results[i].observed - results[i].estimated;
    sse += residual*residual;
  }
  return sse;
}

/// The mean of squared errors.  The response variable at dataSet's 
/// defaultIndex is interpreted as the true function value.
double Surface::mse(SurfData& dataSet)
{
  return sse(dataSet) / static_cast<double>(dataSet.size());
}

/// The maximum relative absolute error is computed by dividing the maximum
/// absolute error by the standard deviation of the data.  The response 
/// variable at dataSet's defaultIndex is interpreted to be the true 
/// function value.
double Surface::mrae(SurfData& dataSet)
{
  vector<surfpack::ErrorStruct> results;
  vector<double> trueVals(dataSet.size());
  getValue(dataSet,results);
  double max = abs(results[0].observed - results[0].estimated);
  for (unsigned i = 1; i < results.size(); i++) {
    double curr = abs(results[i].observed -results[i].estimated);
    trueVals[i] = results[i].observed;
    //cout << "residual: " << residual << " sq: " << residual*residual << endl;
    if (curr > max) max = curr;
  }
  
  return max / surfpack::sample_sd(trueVals);
}

// ____________________________________________________________________________
// Commands 
// ____________________________________________________________________________

/// Associates a data set with this Surface object.  When Surface::build()
/// is invoked, this is the data that will be used to create the Surface
/// approximation.
void Surface::setData(SurfData* sd_)
{
  // If this Surface is already listening to a data set, it needs to notify that
  // SurfData object that it is going to stop listening.
  if (this->sd) {
    this->sd->removeListener(this);
  }
  this->sd = sd_;
  // Now request to be added to the new SurfData object's list of listeners.
  if (this->sd) {
    this->sd->addListener(this);
  }
  if (builtOK) {
    dataModified = true;
  }
  responseIndex = sd ? sd->getDefaultIndex() : 0;
}
  
/// Causes the data to be scaled along each dimension so that all of the
/// values lie on the interval [0,1].  scaled_val = (old_val - min_val) /
/// (max_val - min_val).
void Surface::scaleUniform()
{
  this->scaler = new SurfScaler();
}

/// Causes data not to be scaled at all before building
void Surface::noScale()
{
  delete scaler;
  scaler = 0;
}

/// Set the state of the SurfData object to use the same defaultIndex and 
/// set of excludedPoints that were used when the Surface approximation was
/// built
void Surface::prepareData()
{
  if (sd) {
    sd->setExcludedPoints(excludedPoints);
    sd->setDefaultIndex(responseIndex);
  } 
}

/// Return the data set pointed to by member sd if parameter dataSet is NULL.
/// If dataSet is not NULL, then return the SurfData it points to.  If both
/// the parameter dataSet and member sd are NULL, throw an exception.  This
/// method is primarily used with the fitness metrics that can be set to use
/// by default the data that the approximation was created with but can also
/// use another set of data provided by the client.
SurfData& Surface::checkData(SurfData* dataSet)
{
    if (dataSet) {
      return *dataSet;
    } else if (sd) {
      prepareData();
      return *sd;
    } else {
      ostringstream errormsg;
      errormsg << "In Surface::checkData: No data was passed in "
	       << "and this surface has no data." << endl;
      throw SurfData::bad_surf_data(errormsg.str());
    }
}

/// Invoked by data member sd when the data changes 
void Surface::notify(int msg)
{
  if (msg == SurfData::GOING_OUT_OF_EXISTENCE) {
    sd = 0;
    if (builtOK) {
      dataModified = true;
    }
    responseIndex = sd ? sd->getDefaultIndex() : 0;
  } else if (msg == SurfData::DATA_MODIFIED) {
    if (builtOK) {
      dataModified = true;
    }
  }
}
/// Check to make sure that data are acceptable and then build.
/// Do not build if the surface has already been built and the data have not
/// changed.
void Surface::createModel(SurfData* surfData)
{
  // If they passed in a non-NULL data set, use it to build the Surface.
  // Otherwise, use the data already pointed to by member sd.
  if (surfData) {
    setData(surfData);
  }
  if (builtOK && !dataModified) {
    // Model is already valid and will not be rebuilt because the data has not
    // changed.
    return;
  } 
  acceptableData();
  xsize = sd->xSize();
  // If a SurfScaler object has been created, the data should be scaled
  // before building the surface
  if (scaler) {
    sd->setScaler(scaler);
  }
  SurfData& sdRef = *sd;
  build(sdRef);
  excludedPoints = sd->getExcludedPoints();
  builtOK = true;
  dataModified = false;
  if (scaler) {
    sd->setScaler(0);
  }
}

// ____________________________________________________________________________
// I/O 
// ____________________________________________________________________________

/// Write the surface out to a file.  Files with extension .txt are written
/// out in text mode; .srf files are written in binary format (unformatted). 
void Surface::write(const string filename)
{
  const string nameOfSurface = surfaceName();
  // testFileExtension returns true for .srf, false for .txt
  bool binary = testFileExtension(filename); 
  ofstream outfile(filename.c_str(), 
    (binary ? ios::out|ios::binary : ios::out));
  if (!outfile) {
    throw surfpack::file_open_failure(filename);
  } else if (binary) {
    // write out the surface name
    unsigned nameSize = nameOfSurface.size();
    outfile.write(reinterpret_cast<char*>(&nameSize),sizeof(nameSize));
    outfile.write(nameOfSurface.c_str(),nameSize);
    // write out the surface 'details'
    writeBinary(outfile);
    // If data member sd points to the data used to create this surface, write
    // it out as well.
    if (hasOriginalData()) {
      outfile.write(reinterpret_cast<char*>(&responseIndex),
        sizeof(responseIndex));
    } else {
      // a -1 signifies that there is no data
      int dummy = -1;
      outfile.write(reinterpret_cast<char*>(&dummy), sizeof(dummy));
    }
  } else {
    // write out the surface name
    outfile << nameOfSurface << endl;
    // write out the surface 'details'
    writeText(outfile);
    // If data member sd points to the data used to create this surface, write
    // it out as well.
    if (hasOriginalData()) {
      outfile << responseIndex << " response index for surface data" << endl;
    } else {
      // a -1 signifies that there is no data
      outfile << "-1 data for surface not included" << endl;
    }
  }
  if (hasOriginalData()) {
    prepareData();
    if (binary) {
      sd->writeBinary(outfile);
    } else {
      sd->writeText(outfile);
    }
  }
  outfile.close();
}

/// Read the surface from a file.  Files with extension .txt are read in text 
/// mode; others are read in binary format (unformatted).
void Surface::read(const string filename)
{
  // First, open the file
  // testFileExtension returns true for .srf, false for .txt
  bool binary = testFileExtension(filename);
  ifstream infile(filename.c_str(), (binary ? ios::in|ios::binary : ios::in));
  int index;
  if (!infile) {
    throw surfpack::file_open_failure(filename);
  } 
  // The first thing in the file is the name of the surface.
  string nameInFile = surfpack::readName(infile, binary);
  if (nameInFile != surfaceName() ) {
    ostringstream errormsg;
    errormsg << "Bad surface name in file " << filename << ".  "
             << "Expected: " << surfaceName() << "; found: " 
             <<  surfpack::surfaceName(filename) 
             << ".  Cannot build surface." << endl;
    throw surfpack::io_exception(errormsg.str());
  } else if (binary) {
    // read the surface details
    readBinary(infile);
    infile.read(reinterpret_cast<char*>(&index),sizeof(index));
  } else {
    // read the surface details
    string sline;
    readText(infile);
    getline(infile, sline);
    istringstream streamline(sline);
    streamline >> index;
  }
  if (index >= 0) {
    SurfData* fileData = new SurfData(infile, binary);
    this->setData(fileData);
    responseIndex = static_cast<unsigned>(index);
  }    
  infile.close();
  builtOK = true;
  dataModified = false;
}

/// Return true if filename has .srf extension, false if filename has .txt
/// extension.  If neither, throw surfpack::io_exception.
bool Surface::testFileExtension(const std::string& filename) const
{
  if (surfpack::hasExtension(filename,".srf")) {
    return true;
  } else if (surfpack::hasExtension(filename,".txt")) {
    return false;
  } else {
    throw surfpack::io_exception(
      "Unrecognized filename extension.  Use .srf or .txt"
    );
  }
}

// Write the surface out to a stream in text format
ostream& operator<<(ostream& os,Surface& surface)
{ 
  surface.writeText(os); 
  return os;
}
// ____________________________________________________________________________
// Testing
// ____________________________________________________________________________
