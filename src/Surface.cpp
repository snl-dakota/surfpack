// These lines were automatically prepended
#ifndef SURFPACK_CONFIG_H
#define SURFPACK_CONFIG_H
#include "config.h"
#endif
// End of prepended lines

// ----------------------------------------------------------
// Project: SURFPACK++
//
// File:        Surface.cpp
// Author:      Tony Giunta
// Modified:    Eric Cyr 
// Modified:    Mark Richards
//
// Description: 
// + The Surface class provides an interface for the
//   different surface types, this class cannot be instantiated
// ----------------------------------------------------------

#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <set>
#include <string>
#include "surfpack.h"
#include "Surface.h"
#include "SurfPoint.h"
#include "SurfData.h"
#include <cstdlib>
#include <cmath>
//#include "SurfDataIterator.h"
//#include "SkipSurfDataIterator.h"


using namespace std;

// ____________________________________________________________________________
// Creation, Destruction, Initialization 
// ____________________________________________________________________________

Surface::Surface() : sd(0)
{
#ifdef __TESTING_MODE__
  constructCount++;
#endif
  init();
}

Surface::Surface(SurfData* sd) : sd(sd)
{ 
#ifdef __TESTING_MODE__
  constructCount++;
#endif
  init();
  //if(surfData) {
  //  surfData.addListener(this); 
  //}
}

Surface::Surface(const Surface& other) : xsize(other.xsize), 
  builtOK(other.builtOK), dataAdded(other.dataAdded), 
  dataModified(other.dataModified), originalData(other.originalData),
  responseIndex(other.responseIndex)
{
#ifdef __TESTING_MODE__
  constructCount++;
#endif
  this->sd = other.sd ? new SurfData(*(other.sd)) : 0;
}

/// destructor, just removes itself from any SurfData objects
Surface::~Surface() 
{ 
  //if(surfData) {
  //   surfData.removeListener(this);
  //}
#ifdef __TESTING_MODE__
  destructCount++;
#endif
}

void Surface::init()
{
  dataAdded = dataModified = builtOK = 0;
  xsize = sd ? sd->xSize() : 0;
  responseIndex = 0;
}

// Copy constructor goes here

// ____________________________________________________________________________
// Queries 
// ____________________________________________________________________________

/// Return dimensionality of the surface or zero if not built
unsigned Surface::xSize()
{
  return xsize;
}

bool Surface::hasOriginalData() const
{
  return (sd != 0 && !dataAdded && !dataModified); 
}

bool Surface::acceptableData() const
{
  if (!sd) {
    cerr << "Data unacceptable: there is no data" << endl;
    return false;
  } else if (sd->size() < minPointsRequired()) {
    cerr << "Data unacceptable: there is not enough data " << endl;
    return false;
  }
  return true;
}
  
/// Evaulate the approximation surface at point x and return the value.
/// The point x must have the same dimensionality as this surface's SurfData.
/// Make sure the surface is valid and then call evaluate
double Surface::getValue(const std::vector<double>& x)
{
  if (!builtOK || dataModified) {
    createModel();
  }
  
  // By this point, the model should be built correctly
  if (!builtOK) {
    cerr << "Cannot evaluate point because model was not successfully rebuilt" 
         << endl;
    return 0.0;
  } else {
    return evaluate(x);
  }
  
}

double Surface::getValue(const SurfPoint& sp)
{
  return getValue(sp.X());
}

void Surface::getValue(SurfData& sd, std::vector<ErrorStruct>& pts)
{
  for (unsigned i = 0; i < sd.size(); i++) {
    ErrorStruct es;
    es.estimated = evaluate(sd[i].X());
    es.observed = sd.getResponse(i);
    pts.push_back(es);
  }
}

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

/// Evaluate the empirical model at the points in surfData and output
/// the points and their evaluations to os
//double Surface::test(SurfData& surfData, ostream& os)
//{
//  ofstream outfile("results.txt", ios::out);
//  double error = 0.0;
//  for (unsigned i = 0; i < surfData.size(); i++) {
//    SurfPoint& currentPoint = surfData.Point(i);
//    double estimatedResponse = evaluate(currentPoint.X());
//    
//    // print out the points and deviations to a file
//    currentPoint.writeText(outfile);
//    outfile << setw(20) << estimatedResponse 
//            << setw(20) << (currentPoint.F() - estimatedResponse)
//  	  << endl;	
//    error += (currentPoint.F() - estimatedResponse) * 
//      (currentPoint.F() - estimatedResponse);
//  }
//  outfile.close();
//  return error;
//}    

double Surface::goodnessOfFit(const string metricName, SurfData* surfData)
{
  SurfData& sdRef = checkData(surfData);
  if (metricName == "sse") {
    return sse(sdRef);
  } else if (metricName == "mse") {
    return mse(sdRef);
  } else if (metricName == "mrae") {
    return mse(sdRef);
  } else if (metricName == "rsquared") {
    return rSquared(sdRef);
  } else if (metricName == "press") {
    return press(sdRef);
  } else {
    cout << "No error metric of that type in this class" << endl;
  }
  return 0;
}

double Surface::press(SurfData& dataSet)
{
  /// <= test is used because it must be possible to build the surface
  /// even when one point is removed from dataSet
  if (dataSet.size() <= minPointsRequired()) {
    cerr << "Not enough data to compute PRESS" << endl;
    return 0.0;
  } else {
    // If some of the points in the data set are already being excluded,
    // copy all of the non-excluded data points into a new SurfData
    // object where none of the points are excluded.  This will make
    // it easier to do the "leave one out" process
    bool containsInactives = !dataSet.getExcludedPoints().empty();
    SurfData surfData = (containsInactives) 
      ? dataSet.copyActive() : dataSet;
    double pressValue = 0.0;
    unsigned i = 0;

    // For each data point, make a new surface that has all of the
    // data points except the current one.  Evaluate the new
    // surface at the omitted point and compute the residual
    // between the actual value and the value predicted by the new 
    // model at that point.
    unsigned totalPoints = surfData.size();
    i = 0;
    set<unsigned> pointToSkip;
    while (i < totalPoints) {
      pointToSkip.clear();
      surfData.setExcludedPoints(pointToSkip);
      // Get the point from the full data set corresponding to the index that 
      // will be skipped.  (This point will not be available to test after it
      // has been marked for exclusion
      const SurfPoint& currentPoint = surfData[i];
      pointToSkip.insert(i);
      surfData.setExcludedPoints(pointToSkip);

      Surface* allButOne = makeSimilarWithNewData(&surfData);
      double fTilde = allButOne->getValue(currentPoint.X());
      // actual value of the function at this point
      double fx = currentPoint.F();
      // add the square of the residual (fTilde - fx)
      pressValue += (fTilde - fx) * (fTilde - fx);
      delete allButOne;
      allButOne = 0;
      i++;
      //cout << setw(5) << setprecision(0) 
      //     << ( 100.0*static_cast<double>(i)/static_cast<double>(surfData.size()+1) ) 
      //     << "%\r" << flush;
    }
    cout << setprecision(6) << endl;
    pressValue = sqrt(pressValue/static_cast<double>(surfData.size()+1));
    return pressValue;
  }
  return 0.0;
}

double Surface::rSquared(SurfData& dataSet)
{
  double sumObserved = 0.0;
  double sumOfSquaresObserved = 0.0;
  double residualSumOfSquares = 0.0;
  double totalSumOfSquares = 0.0;
  for (unsigned i = 0; i < dataSet.size(); i++) {
    double observedF = dataSet.getResponse(i);
    double estimatedF = evaluate(dataSet[i].X());
    double residual = observedF - estimatedF;
    //cout << setw(4) << i
    //     << setw(15) << observedF
    //     << setw(15) << estimatedF
    //     << setw(15) << residual
    //     << setw(15) << residual*residual
    //     << endl;
    residualSumOfSquares += residual * residual;
    sumObserved += observedF;
    sumOfSquaresObserved += observedF * observedF;
  }
  totalSumOfSquares = sumOfSquaresObserved - 
    (sumObserved * sumObserved / dataSet.size());
  double rSquaredValue = 1.0 - residualSumOfSquares / totalSumOfSquares;
  //cout << "totalSumOfSquares: " << totalSumOfSquares << endl;
  //cout << "sumOfSquaresObserved: " << sumOfSquaresObserved << endl;
  //cout << "sumObserved: " << sumObserved << " " << sumObserved*sumObserved << endl;
  //cout << "count: " << dataSet.size() << endl;
  //cout << "residualSumOfSquares: " << residualSumOfSquares << endl;
  //cout << "rSquared: " << rSquaredValue << endl;
                                                                                                                 

  return (rSquaredValue < 0) ? 0 : rSquaredValue;
} 
      
double Surface::sse(SurfData& dataSet)
{
  vector<ErrorStruct> results;
  getValue(dataSet,results);
  double sse = 0.0;
  for (unsigned i = 0; i < results.size(); i++) {
    double residual = results[i].observed - results[i].estimated;
    //cout << "residual: " << residual << " sq: " << residual*residual << endl;
    sse += residual*residual;
  }
  return sse;
}

double Surface::mrae(SurfData& dataSet)
{
  vector<ErrorStruct> results;
  vector<double> trueVals(dataSet.size());
  getValue(dataSet,results);
  double max = abs(results[0].observed - results[0].estimated);
  for (unsigned i = 1; i < results.size(); i++) {
    double curr = abs(results[i].observed -results[i].estimated);
    trueVals[i] = results[i].observed;
    //cout << "residual: " << residual << " sq: " << residual*residual << endl;
    if (curr > max) max = curr;
  }
  
  return max / sample_sd(trueVals);
}

double Surface::mse(SurfData& dataSet)
{
  return sse(dataSet) / static_cast<double>(dataSet.size());
}

//double Surface::computePressStatistic()
//{
//  double fpress = 0;
//  //DataIterator dataIterator;
//  //surfData.getIterator(dataIterator);
//  //const SurfPoint* currentPoint = 0;
//  unsigned i = 0;
//  SkipSurfDataIterator skipIt(surfData);
//  while (i < surfData.size()) {
//    // leave out one point from the data set
//    const SurfPoint currentPoint = surfData.getPoint(i); 
//    assert(currentPoint);
//    skipIt.unSkipAll();
//    skipIt.skipPoint(i);
//    //cout << "Skipping: " << i << endl;
//    build(&skipIt);
//    double fTilde = evaluate(currentPoint.getX());
//    // actual value of the function at this point
//    double fx = currentPoint.getF();
//    // add the square of the deviation fTilde - fx
//    fpress += (fTilde - fx) * (fTilde - fx);
//    //dataIterator.nextElement();
//    i++;
//  }
//  skipIt.unSkipAll();
//  build(&skipIt);
//  fpress = sqrt(fpress/i);
//  return fpress;
//}

// ____________________________________________________________________________
// Commands 
// ____________________________________________________________________________

/// Associates a data set with this surface object.  If this surface has
/// already been built, it is invalidated
void Surface::setData(SurfData* sd)
{
  //Do a sanity check on the data
  this->sd = sd;
  dataModified = true;
  xsize = sd ? sd->xSize() : 0;
}
  
/// Set the state of the SurfData object to use the default index and points
/// associated with this surface
void Surface::prepareData()
{
  if (sd) {
    sd->setExcludedPoints(excludedPoints);
    sd->setDefaultIndex(responseIndex);
  } else {
    cerr << "There is no data to prepare." << endl;
  }
}

/// Checks to make sure the data passed in is not null.  If it is, sets it 
/// to point to the SurfData used to create the object.  If that is also
/// non-existent, it is an error.  
SurfData& Surface::checkData(SurfData* dataSet)
{
  if (!dataSet) {
    if (sd) {
      prepareData();
      return *sd;
    } else {
      cerr << "In Surface::checkData: No data was passed in and this surface has"
           << " no data." << endl;
    }
  }
  return *dataSet;
}

/// Check to make sure that data are acceptable and then build.
/// Do not build if the surface has already been built and the data have not
/// changed
void Surface::createModel(SurfData* surfData)
{
  if (surfData) {
    setData(surfData);
  }
  if (builtOK && !dataAdded && !dataModified) {
    cerr << "Model is already valid and will not be rebuilt because the"
         << " data set has not changed"
         << endl;
    return;
  } 

  if (acceptableData()) {
    SurfData& sdRef = *sd;
    build(sdRef);
    excludedPoints = sd->getExcludedPoints();
    responseIndex = sd->getDefaultIndex();
    originalData = builtOK = true;
    dataAdded = dataModified = false;
  } else {
    cerr << "Model could not be built because the data was not acceptable" 
         << endl;
  }
}

// ____________________________________________________________________________
// I/O 
// ____________________________________________________________________________

void Surface::write(const string filename)
{
  const string nameOfSurface = surfaceName();
  bool binary = (filename.find(".txt") != filename.size() - 4);
  ofstream outfile(filename.c_str(), (binary ? ios::out|ios::binary : ios::out));
  if (!outfile) {
    cerr << "File named \"" 
         << filename
	 << "\" could not be opened." 
	 << endl;
    return;
  } else if (binary) {
    // write out the surface name
    unsigned nameSize = nameOfSurface.size();
    outfile.write(reinterpret_cast<char*>(&nameSize),sizeof(nameSize));
    outfile.write(nameOfSurface.c_str(),nameSize);
    // write out the surface 'details'
    writeBinary(outfile);
    if (sd) {
      outfile.write(reinterpret_cast<char*>(&responseIndex),
        sizeof(responseIndex));
    } else {
      int dummy = -1;
      outfile.write(reinterpret_cast<char*>(&dummy), sizeof(dummy));
    }
  } else {
    // write out the surface name
    outfile << nameOfSurface << endl;
    // write out the surface 'details'
    writeText(outfile);
    if (sd) {
      outfile << responseIndex << " response index for surface data" << endl;
    } else {
      outfile << "-1 data for surface not included" << endl;
    }
  }
  if (originalData) {
    prepareData();
    if (binary) {
      sd->writeBinary(outfile);
    } else {
      sd->writeText(outfile);
    }
  }
  outfile.close();
}

void Surface::read(const string filename)
{
  bool binary = (filename.find(".txt") != filename.size() - 4);
  ifstream infile(filename.c_str(), (binary ? ios::in|ios::binary : ios::in));
  int index;
  if (!infile) {
    cerr << "File named \"" 
         << filename
	 << "\" could not be opened." 
	 << endl;
    return;
  } else if (binary) {
    // read surface name
    unsigned nameSize;
    infile.read(reinterpret_cast<char*>(&nameSize),sizeof(nameSize));
    char* surfaceType = new char[nameSize+1];
    infile.read(surfaceType,nameSize);
    surfaceType[nameSize] = '\0';
    string nameInFile(surfaceType);
    delete [] surfaceType;
    if (nameInFile != surfaceName()) {
      cerr << "Surface name in file is not " << surfaceName() << "." << endl;
      cerr << "Cannot build surface." << endl;
      return;
    }
    // read the surface details
    readBinary(infile);
    infile.read(reinterpret_cast<char*>(&index),sizeof(index));
  } else {
    // read surface name 
    string sline;
    getline(infile,sline);
    if (sline!= surfaceName()) {
      cerr << "Surface name in file is not " << surfaceName() << "." << endl;
      cerr << "Cannot build surface." << endl;
      return;
    }
    // read the surface details
    readText(infile);
    getline(infile, sline);
    istringstream streamline(sline);
    streamline >> index;
  }
  if (index >= 0) {
    sd = new SurfData(infile, binary);
    responseIndex = static_cast<unsigned>(index);
    originalData = true;
  }    
  infile.close();
  builtOK = true;
  dataAdded = dataModified = false;
}

ostream& operator<<(ostream& os,Surface& surface)
{ 
  surface.writeText(os); 
  return os;
}


/// Write the associated data to a stream.  Not all iterators will use all of 
/// the data available in their SurfData objects, so the writing must 
/// necessarily go through the iterator.
//void Surface::writeData(std::ostream& os, bool binary)
//{
//  if (!dataItr) {
//    cerr << "No data to write in Surface::writeData" << endl;
//  } else {
//    unsigned s = dataItr->elementCount();
//    unsigned xsize = dataItr->xSize();
//    unsigned fsize = dataItr->currentElement().fSize();
//    // write out header (number of points, dimension, number of responses)
//    if (binary) {
//      os.write((char*)&s,sizeof(s));
//      os.write((char*)&xsize,sizeof(xsize));
//      os.write((char*)&fsize,sizeof(fsize));
//    } else {
//      os << s << " data points" << endl
//         << xsize << " input variables" << endl
//       	 << fsize << " response variables" << endl;
//    }
//    // write out each point, one at a time
//    dataItr->toFront();
//    for (unsigned i = 0; i < dataItr->elementCount(); i++) {
//      if (binary) {
//	dataItr->currentElement().writeBinary(os);
//      } else {
// 	dataItr->currentElement().writeText(os);
//      }
//      dataItr->nextElement();
//    }
//  }
//}

/// Read the data from a file and create a SurfDataIterator wrapper for them
//void Surface::readData(std::istream& is, bool binary)
//{
//  delete dataItr;
//  sd = new SurfData(is, binary);
//  dataItr = new SurfDataIterator(*sd);
//}

// ____________________________________________________________________________
// Testing
// ____________________________________________________________________________


#ifdef __TESTING_MODE__
  int Surface::constructCount = 0;
  int Surface::destructCount = 0;
#endif

