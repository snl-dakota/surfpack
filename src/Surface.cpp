/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#include "surfpack_config.h"

#include "surfpack.h"
#include "SurfData.h"
#include "Surface.h"
#include "SurfScaler.h"

using namespace std;
using surfpack::dbg;
const int dbgsrf = 0;

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
  init();
  setData(sd_);
}

/// Deep copies are made of most data members.  For the sd member, the new
/// object simply asks to be added to sd's list of listeners. 
// sd must be initialized to NULL before setData(...) is called, because 
// setData(...) will call delete on the old value of sd, which will be
// garbage if it's not initialized to NULL.
Surface::Surface(const Surface& other) 
  : sd(0),  scaler(0), xsize(other.xsize), 
  builtOK(other.builtOK), dataModified(other.dataModified), 
  responseIndex(other.responseIndex)
{
  if (other.scaler) {
    scaler = new SurfScaler(*other.scaler);
  }
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
unsigned Surface::xSize() const
{
  return this->xsize;
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
      errormsg << "Data unacceptable: this surface requires "
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
double Surface::getValue(const vector<double>& x)
{
  if (!builtOK || dataModified) {
    createModel();
  }
  if (!scaler) {
    return evaluate(x);
  } else {
    SurfPoint sp(x);
    sp.setScaler(scaler);
    return scaler->descaleResponse(responseIndex,
      evaluate(sp.X()));
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
void Surface::getValue(SurfData& surf_data, vector<surfpack::ErrorStruct>& pts)
{
  for (unsigned i = 0; i < surf_data.size(); i++) {
    surfpack::ErrorStruct es;
    es.estimated = getValue(surf_data[i].X());
    es.observed = surf_data.getResponse(i);
    pts.push_back(es);
  }
}

/// Evaluate the approximation surface at each point in the parameter
/// SurfData object. Return lists of the observed values (retrieved via
/// sd.getResponse() for each point sd) and corresponding predicted values 
/// returned by sd.getResponse()) and the estimated value.  These lists are
/// returned through the second and third parameters.
void Surface::getValue(SurfData& surf_data, vector<double>& observed_vals,
  vector<double>& predicted_vals)
{
  observed_vals.resize(surf_data.size());
  predicted_vals.resize(surf_data.size());
  for (unsigned i = 0; i < surf_data.size(); i++) {
    predicted_vals[i]  = getValue(surf_data[i].X());
    observed_vals[i] = surf_data.getResponse(i);
  }
}

/// Evaluate the approximation surface at each point in the parameter
/// surfData object.  Append the evaluations as a new response variable in
/// the data set.
unsigned Surface::getValue(SurfData& surfData)
{
  set<unsigned> emptySet;
  surfData.setExcludedPoints(emptySet);
  vector<double> newValues;
  for (unsigned i = 0; i < surfData.size(); i++) {
    newValues.push_back(getValue(surfData[i].X()));
  }
  return surfData.addResponse(newValues);
}

/// Modify one of the configuration options for this surface type 
void Surface::config(const Arg& arg)
{
  // Other arguments are processed by config methods in the child class.
  // If the child class sees an argument it does not recognize, it should
  // call this method.
  if (arg.name == "response_index") {
    unsigned index = static_cast<unsigned>(arg.getRVal()->getInteger());
    if (!sd) {
      throw string("Cannot set response index on NULL data");
    } else {
      sd->setDefaultIndex(index);
      this->responseIndex = index; 
    }
  } else if (arg.name == "scaling") {
    string scaleOption = arg.getRVal()->getIdentifier();
    if (scaleOption == "uniform") {
      scaleUniform();
    } else if (scaleOption == "auto") {
      //autoScale();
    } else if (scaleOption == "none") {
      noScale();
    } else {
      throw string("Unrecognized option for surface parameter 'scaling'");
    }  
  } else if (arg.name == "norm_scale" || arg.name == "log_scale") {
    scalingArg(arg);
  } else if (arg.name == "xsize") {
    setXSize(arg.getRVal()->getInteger());
  }
}

void Surface::gradient(const vector<double> & x, 
vector<double>& gradient_vector)
{
  throw string("This surface type does not support hessians");
}

void Surface::hessian(const vector<double> & x, 
SurfpackMatrix<double>& hessian)
{
  throw string("This surface type does not support hessians");
}

/// Process a list of configuration options
void Surface::configList(const ArgList& arglist)
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
  // Make sure the surfData passed in is not null, and that it's configured
  // to use the correct points and response value for this surface
  SurfData& sdRef = checkData(surfData);
  if (metricName == "rsquared") {
    return rSquared(sdRef);
  } else if (metricName == "press") {
    return nFoldCrossValidation(sdRef,sdRef.size());
  } else {
    // The rest of these metrics all have many computations in common
    // and are grouped together in the genericMetric method
    vector<double> observed;
    vector<double> predicted;
    getValue(sdRef, observed, predicted);
    if (metricName == "min_abs" ) {
      return genericMetric(observed,predicted,MT_MINIMUM,ABSOLUTE);
    } else if (metricName == "max_abs") {
      return genericMetric(observed,predicted,MT_MAXIMUM,ABSOLUTE);
    } else if (metricName == "sum_abs") {
      return genericMetric(observed,predicted,MT_SUM,ABSOLUTE);
    } else if (metricName == "mean_abs") {
      return genericMetric(observed,predicted,MT_MEAN,ABSOLUTE);
    } else if (metricName == "max_relative") {
      return genericMetric(observed,predicted,MT_RELATIVE_MAXIMUM,ABSOLUTE);
    } else if (metricName == "mean_relative") {
      return genericMetric(observed,predicted,MT_RELATIVE_AVERAGE,ABSOLUTE);
    } else if (metricName == "min_squared" ) {
      return genericMetric(observed,predicted,MT_MINIMUM,SQUARED);
    } else if (metricName == "max_squared") {
      return genericMetric(observed,predicted,MT_MAXIMUM,SQUARED);
    } else if (metricName == "sum_squared") {
      return genericMetric(observed,predicted,MT_SUM,SQUARED);
    } else if (metricName == "mean_squared") {
      return genericMetric(observed,predicted,MT_MEAN,SQUARED);
    } else if (metricName == "min_scaled" ) {
      return genericMetric(observed,predicted,MT_MINIMUM,SCALED);
    } else if (metricName == "max_scaled") {
      return genericMetric(observed,predicted,MT_MAXIMUM,SCALED);
    } else if (metricName == "sum_scaled") {
      return genericMetric(observed,predicted,MT_SUM,SCALED);
    } else if (metricName == "mean_scaled") {
      return genericMetric(observed,predicted,MT_MEAN,SCALED);
    } else if (metricName == "root_mean_squared") {
      return rootMeanSquared(observed,predicted);
    } else {
      throw string("No error metric of that type in this class");
    }
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
    // model at that point.
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

unsigned block_low(unsigned i, unsigned n, unsigned p) { return i*n/p;}
unsigned block_high(unsigned i, unsigned n, unsigned p) { return (i+1)*n/p-1;}
double Surface::nFoldCrossValidation(SurfData& data, unsigned n)
{
  double total_error = 0.0;
  // If n is evenly divisible by data.size(), then the number of points that
  // gets left out each time is the same: n/data.size().  We need to make sure that
  // even when the maximum number of points is excluded, there are still enough
  // points to run the algorithm
  unsigned data_size = data.size();
  unsigned max_exclusions = data_size / n; 
  if (data_size % n != 0) max_exclusions++;
  if ((data_size - max_exclusions) < minPointsRequired()) {
    throw SurfData::
      bad_surf_data("Not enough data to compute n-fold cross validation.");
  } else {
    // If some of the points in the data set are already being excluded,
    // copy all of the non-excluded data points into a new SurfData
    // object where none of the points are excluded.  This will make
    // the rest of the process simpler and more intuitive, albeit slightly less
    // efficient 
    bool containsInactives = !data.getExcludedPoints().empty();
    SurfData active_set = (containsInactives) 
      ? data.copyActive() : data;
    data_size = active_set.size();

    // Now create the partitions
    vector<unsigned> skip_points(active_set.size());
    for (unsigned i = 0; i < active_set.size(); i++) skip_points[i] = i;
    ///\todo manage the seed and random number generation
    surfpack::shared_rng().seed(0);
    random_shuffle(skip_points.begin(),skip_points.end(),surfpack::shared_rng());

    if (dbgsrf) copy(skip_points.begin(),skip_points.end(),
      ostream_iterator<unsigned>(cout," "));
    // Compute the error for each partition
    for (unsigned part_index = 0; part_index < n; part_index++) {
      // Determine which points to exclude for the current partition
      set<unsigned> points_to_exclude;
      unsigned lower_index = block_low(part_index,data_size,n);
      unsigned upper_index = block_high(part_index,data_size,n);
      dbg(dbgsrf) << "Partition " << part_index << " ex. pts.:" ;
      for (unsigned k = lower_index; k <= upper_index; k++) {
        points_to_exclude.insert(skip_points[k]);
        dbg(dbgsrf) << skip_points[k] << " ";
      }
      dbg(dbgsrf) << '\n';
      // debug code
      active_set.setExcludedPoints(points_to_exclude);
      dbg(dbgsrf) << "active_set.size(): " << active_set.size() << "\n";
      Surface* current_surf = makeSimilarWithNewData(&active_set);
      current_surf->createModel();
      double partition_error = 0.0;
      // Now evaluate the surface at the points with known values that were 
      // excluded. To make sure that the indices of the points are coorespond
      // to the right points, we need to unexclude all points.
      
      // Accumulate the squared residuals
      set<unsigned> no_exclusions;
      active_set.setExcludedPoints(no_exclusions);
      for (unsigned k = lower_index; k <= upper_index; k++) {
        double observed = active_set[skip_points[k]].F(this->responseIndex);
        double predicted = 
          current_surf->getValue(active_set[skip_points[k]].X());
        dbg(dbgsrf) << "resid " << skip_points[k] << ": " << observed << " " << predicted 
	     << '\n';
        partition_error += (observed - predicted)*(observed-predicted);
      }
      total_error += partition_error;
      delete current_surf;
    }
    // To get a value analogous to the root-mean-squared value (instead of
    // something analogous to SSE, uncomment the following line
    //total_error = sqrt(total_error/active_set.size());
  }
  return total_error;
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
      
/// Compute one of several goodness of fit metrics.  The observed parameter
/// should be a list of observed (or true) function values; the vector of
/// predicted values gives the corresponding estimates from this surface.
/// The dt parameter specifies the kind of residuals to compute.  ABSOLUTE
/// residuals are (observed - predicted), SQUARED residuals are the squares 
/// of the absolute residuals.  SCALED residuals are the ABSOLUTE residuals
/// divided by the observed value.  Given the type of residuals, the client
/// may request the min, max, sum, or mean of the set of residuals over all
/// the given data points.  Two additional metrics are possible.  The
/// relative maximum absolute error is the maximum absolute error divided
/// by the standard deviation of the observed data.  The relative average
/// absolute error is the mean absolute error divided by the standard 
/// deviation of observed data. 
double Surface::genericMetric(vector<double>& observed,
  vector<double>& predicted, enum MetricType mt, enum DifferenceType dt)
{
  /// Create a vector for storing the differences (absolute, relative, or 
  /// scaled) between observed and predicted values.
  vector<double> diffs;
  // Iterator for capturing the max or min using STL generic algorithm
  vector<double>::iterator iter;
  // Compute the desired residuals (obs_i - pred_i); they may be
  // absolute, squared, or scaled (obs_i - pred_i)/obs_i, depending
  // on the value of dt
  surfpack::differences(diffs, observed, predicted, dt);
  // The main values for mt are min, max, sum, and mean
  // Two special values give additional metrics: relative maximum and
  // relative average
  switch (mt) {
    // Relative maximum absolute error = max absolute error divided by the
    // standard deviation of the obs_i values.
    case MT_RELATIVE_MAXIMUM:
      iter = max_element(diffs.begin(),diffs.end());
      return *iter / surfpack::sample_sd(observed);
    // Relative average absolute error is the average absolute
    // error (sum_over(obs_i - pred_i) / n) divided by the standard
    // deviation of the obs_i values
    case MT_RELATIVE_AVERAGE:
      return surfpack::sum_vector(diffs) / 
        (observed.size() * surfpack::sample_sd(observed));
    case MT_MINIMUM:
      iter = min_element(diffs.begin(),diffs.end());
      return *iter;
    case MT_MAXIMUM:
      iter = max_element(diffs.begin(),diffs.end());
      return *iter;
    case MT_SUM:
      return surfpack::sum_vector(diffs);
    case MT_MEAN:
    default:
      return surfpack::mean(diffs);
  }
}

double Surface::rootMeanSquared(vector<double>& observed,
  vector<double>& predicted)
{
  return sqrt(genericMetric(observed,predicted,MT_MEAN,SQUARED));
}
// ____________________________________________________________________________
// Commands 
// ____________________________________________________________________________

/// Associates a data set with this Surface object.  When Surface::build()
/// is invoked, this is the data that will be used to create the Surface
/// approximation.
void Surface::setData(SurfData* sd_in)
{
  if (sd_in) {
    dbg(dbgsrf) << "sd_in->size(): " << sd_in->size() << "\n";
  }
    
  // If this Surface is already listening to a data set, it needs to notify that
  // SurfData object that it is going to stop listening.
  if (this->sd) {
    this->sd->removeListener(this);
  }
  this->sd = sd_in;
  // Now request to be added to the new SurfData object's list of listeners.
  if (this->sd) {
    this->sd->addListener(this);
  }
  if (builtOK) {
    dataModified = true;
  }
  responseIndex = sd ? sd->getDefaultIndex() : 0;
  if (sd) xsize = sd->xSize();
}
  
/// Causes the data to be scaled along each dimension so that all of the
/// values lie on the interval [0,1].  scaled_val = (old_val - min_val) /
/// (max_val - min_val).
void Surface::scaleUniform()
{
  assert(sd);
  delete scaler;
  this->scaler = new SurfScaler();
  scaler->normalizeAll(*sd);
}

/// Causes data not to be scaled at all before building
void Surface::noScale()
{
  delete scaler;
  scaler = 0;
}

void Surface::scalingArg(const Arg& arg) 
{
  assert(sd);
  if (!scaler) {
    scaler = new SurfScaler;
  }
  scaler->config(*sd,arg);    
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
  xsize = sd->xSize();
  acceptableData();
  // If a SurfScaler object has been created, the data should be scaled
  // before building the surface
  if (scaler) {
    sd->setScaler(scaler);
    //cout << scaler->asString() << endl;
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

/// Write the surface out to a file.  Files with extension .sps are written
/// out in text mode; .bsps files are written in binary format (unformatted). 
void Surface::write(const string filename)
{
  const string nameOfSurface = surfaceName();
  // hasBinaryFileExtension returns true for .bsps, false for .sps
  bool binary = hasBinaryFileExtension(filename); 
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

/// Read the surface from a file.  Files with extension .sps are read in text 
/// mode; .bsps are read in binary format (unformatted).
void Surface::read(const string filename)
{
  // First, open the file
  // hasBinaryFileExtension returns true for .bsps, false for .sps
  bool binary = hasBinaryFileExtension(filename);
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

/// Return true if filename has .bsps extension, false if filename has .sps
/// extension.  If neither, throw surfpack::io_exception.
bool Surface::hasBinaryFileExtension(const string& filename) const
{
  if (surfpack::hasExtension(filename,".bsps")) {
    return true;
  } else if (surfpack::hasExtension(filename,".sps")) {
    return false;
  } else {
    throw surfpack::io_exception(
      "Unrecognized filename extension.  Use .sps or .bsps"
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
