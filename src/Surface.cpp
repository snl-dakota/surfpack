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
#include <iomanip>
#include <vector>
#include <string>
#include "surfpack.h"
#include "Surface.h"
#include "SurfPoint.h"
#include "SurfData.h"
#include "SurfDataIterator.h"
#include "SkipSurfDataIterator.h"


using namespace std;

// ____________________________________________________________________________
// Creation, Destruction, Initialization 
// ____________________________________________________________________________

Surface::Surface() : valid(false), originalData(false), dataItr(0), sd(0)
{
#ifdef __TESTING_MODE__
  constructCount++;
#endif
}

Surface::Surface(AbstractSurfDataIterator* dataItr) :
  valid(false), originalData(false), dataItr(dataItr), sd(0)
{ 
#ifdef __TESTING_MODE__
  constructCount++;
#endif
  //if(surfData) {
  //  surfData.addListener(this); 
  //}
}

/// destructor, just removes itself from any SurfData objects
Surface::~Surface() 
{ 
  delete dataItr;
  dataItr = 0;
  delete sd;
  sd = 0;
  //if(surfData) {
  //   surfData.removeListener(this);
  //}
#ifdef __TESTING_MODE__
  destructCount++;
#endif
}

// Copy constructor goes here

// ____________________________________________________________________________
// Queries 
// ____________________________________________________________________________

bool Surface::isValid() const
{
  return valid;
}

bool Surface::hasOriginalData() const
{
  return (dataItr != 0); 
}

bool Surface::acceptableData() const
{
  if (!dataItr) {
    return false;
  } else if (dataItr->elementCount() < minPointsRequired()) {
    return false;
  }
  return true;
}

double Surface::evaluate(const SurfPoint& sp)
{
  return evaluate(sp.X());
}

void Surface::evaluate(SurfData& surfData)
{
  unsigned newindex = surfData.addResponse();
  SurfDataIterator surfIt(surfData);
  surfIt.toFront();
  while (!surfIt.isEnd()) {
    SurfPoint& current = surfIt.currentElement();
    double estimatedResponse = evaluate(current.X());
    current.F(newindex, estimatedResponse);
    surfIt.nextElement();
  }
}

/// Evaluate the empirical model at the points in surfData and output
/// the points and their evaluations to os
double Surface::test(SurfData& surfData, ostream& os)
{
  ofstream outfile("results.txt", ios::out);
  double error = 0.0;
  for (unsigned i = 0; i < surfData.size(); i++) {
    SurfPoint& currentPoint = surfData.Point(i);
    double estimatedResponse = evaluate(currentPoint.X());
    
    // print out the points and deviations to a file
    currentPoint.writeText(outfile);
    outfile << setw(20) << estimatedResponse 
            << setw(20) << (currentPoint.F() - estimatedResponse)
  	  << endl;	
    error += (currentPoint.F() - estimatedResponse) * 
      (currentPoint.F() - estimatedResponse);
  }
  outfile.close();
  return error;
}    

double Surface::errorMetric(const string metricName, AbstractSurfDataIterator* itr)
{
  if (metricName == "press") {
    //cout << "Call Press from Surface" << endl;
    return press(itr);
  } else if (metricName == "rsquared") {
    return rSquared(itr);
  } else if (metricName == "sse") {
    return sse(itr);
  } else if (metricName == "mse") {
    return mse(itr);
  } else {
    cout << "No error metric of that type in this class" << endl;
  }
  return 0;
}

double Surface::press(AbstractSurfDataIterator* itr)
{
  ensureValidity();
  if (!itr) {
    if (dataItr) {
      itr = dataItr;
    } else {
      cerr << "Cannot compute PRESS without data" << endl;
      return 0.0;
    }
  }
  if (itr->elementCount() <= 1) {
    cerr << "Not enough data to compute PRESS" << endl;
    return 0.0;
  } else {
    double pressValue = 0.0;
    unsigned i = 0;
    // create SurfData objects with all the data points that were
    // used to create this surface
    vector<SurfPoint> sps;
    itr->toFront();
    for (i = 0; i < itr->elementCount(); i++) {
      sps.push_back(itr->currentElement());
      itr->nextElement();
    }
    SurfData sd(sps);

    // For each data point, make a new surface that has all of the
    // data points except the current one.  Evaluate the new
    // surface at the omitted point and compute the residual
    // between the actual value and the value predicted by the new 
    // model at that point.
    i = 0;
    while (i < sd.size()) {
      // Create new iterator.  Iterator will be deleted by the surface
      // that it is used to build.
      SkipSurfDataIterator* skipIt = 
        new SkipSurfDataIterator(sd, itr->responseIndex());
      // leave out one point from the data set
      skipIt->skipPoint(i);
      const SurfPoint& currentPoint = sd.Point(i);
      //
      Surface* allButOne = makeSimilarWithNewData(skipIt);
      double fTilde = allButOne->evaluate(currentPoint.X());
      // actual value of the function at this point
      double fx = currentPoint.F(itr->responseIndex());
      // add the square of the residual (fTilde - fx)
      pressValue += (fTilde - fx) * (fTilde - fx);
      delete allButOne;
      allButOne = 0;
      i++;
      cout << setw(5) << setprecision(0) 
           << ( 100.0*static_cast<double>(i)/static_cast<double>(sd.size()) ) 
           << "%\r" << flush;
    }
    cout << setprecision(6) << endl;
    pressValue = sqrt(pressValue/static_cast<double>(sd.size()));
    return pressValue;
  }
}

double Surface::rSquared(AbstractSurfDataIterator* itr)
{
  ensureValidity();
  double rSquaredValue = 0.0;
  if (!itr) {
    if (dataItr) {
      itr = dataItr;
    } else {
      cerr << "Cannot compute rSquared without data" << endl;
      return 0.0;
    }
  }
  double sumObserved = 0.0;
  double sumOfSquaresObserved = 0.0;
  double residualSumOfSquares = 0.0;
  double totalSumOfSquares = 0.0;
  itr->toFront();
  for (unsigned i = 0; i < itr->elementCount(); i++) {
    SurfPoint& sp = itr->currentElement();
    double observedF = sp.F(itr->responseIndex());
    double estimatedF = evaluate(sp.X());
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
    itr->nextElement();
  }
  totalSumOfSquares = sumOfSquaresObserved - 
    (sumObserved * sumObserved / itr->elementCount());
  rSquaredValue = 1.0 - residualSumOfSquares / totalSumOfSquares;
  //cout << "totalSumOfSquares: " << totalSumOfSquares << endl;
  //cout << "sumOfSquaresObserved: " << sumOfSquaresObserved << endl;
  //cout << "sumObserved: " << sumObserved << " " << sumObserved*sumObserved << endl;
  //cout << "count: " << itr->elementCount() << endl;
  //cout << "residualSumOfSquares: " << residualSumOfSquares << endl;
  //cout << "rSquared: " << rSquaredValue << endl;
  
  return (rSquaredValue < 0) ? 0 : rSquaredValue;
} 
      
double Surface::sse(AbstractSurfDataIterator* itr)
{
  ensureValidity();
  if (!itr) {
    if (dataItr) {
      itr = dataItr;
    } else {
      cerr << "Cannot compute sse without data" << endl;
      return 0.0;
    }
  }
  vector<ErrorStruct> results;
  evaluate(itr,results);
  double sse = 0.0;
  for (unsigned i = 0; i < results.size(); i++) {
    double residual = results[i].observed - results[i].estimated;
    cout << "residual: " << residual << " sq: " << residual*residual << endl;
    sse += residual*residual;
  }
  return sse;
}

double Surface::mse(AbstractSurfDataIterator* itr)
{
  if (!itr) {
    if (dataItr) {
      itr = dataItr;
    } else {
      cerr << "Cannot compute sse without data" << endl;
      return 0.0;
    }
  }
  return sse(itr) / static_cast<double>(itr->elementCount());
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

void Surface::ensureValidity() 
{
  if(!valid) {
      build();
  }
}

void Surface::evaluate(AbstractSurfDataIterator* itr, 
    std::vector<ErrorStruct>& pts)
{
  if (!itr) {
    cerr << "No data in Surface::evaluate" << endl;
    return;
  }
  itr->toFront();
  for (unsigned i = 0; i < itr->elementCount(); i++) {
    SurfPoint& sp = itr->currentElement();
    ErrorStruct es;
    es.estimated = evaluate(sp.X());
    es.observed = sp.F(itr->responseIndex());
    pts.push_back(es);
    itr->nextElement();
  }
}

// ____________________________________________________________________________
// I/O 
// ____________________________________________________________________________

void Surface::write(const string filename)
{
  bool binary = (filename.find(".txt") != filename.size() - 4);
  ofstream outfile(filename.c_str(), (binary ? ios::out|ios::binary : ios::out));
  const string nameOfSurface = surfaceName();
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
  } else {
    // write out the surface name
    outfile << nameOfSurface << endl;
    // write out the surface 'details'
    writeText(outfile);
  }
  writeData(outfile, binary);
  outfile.close();
}

void Surface::read(const string filename)
{
  bool binary = (filename.find(".txt") != filename.size() - 4);
  ifstream infile(filename.c_str(), (binary ? ios::in|ios::binary : ios::in));
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
      cerr << "Surface name in file is not 'Polynomial'." << endl;
      cerr << "Cannot build surface." << endl;
      return;
    }
    // read the surface details
    readBinary(infile);
  } else {
    // read surface name 
    string nameInFile;
    getline(infile,nameInFile);
    if (nameInFile != surfaceName()) {
      cerr << "Surface name in file is not 'Polynomial'." << endl;
      cerr << "Cannot build surface." << endl;
      return;
    }
    // read the surface details
    readText(infile);
  }
  readData(infile, binary);
  originalData = true;
  infile.close();
}

ostream& operator<<(ostream& os,Surface& surface)
{ 
  surface.writeText(os); 
  return os;
}


/// Write the associated data to a stream.  Not all iterators will use all of 
/// the data available in their SurfData objects, so the writing must 
/// necessarily go through the iterator.
void Surface::writeData(std::ostream& os, bool binary)
{
  if (!dataItr) {
    cerr << "No data to write in Surface::writeData" << endl;
  } else {
    unsigned s = dataItr->elementCount();
    unsigned xsize = dataItr->xSize();
    unsigned fsize = dataItr->currentElement().fSize();
    // write out header (number of points, dimension, number of responses)
    if (binary) {
      os.write((char*)&s,sizeof(s));
      os.write((char*)&xsize,sizeof(xsize));
      os.write((char*)&fsize,sizeof(fsize));
    } else {
      os << s << " data points" << endl
         << xsize << " input variables" << endl
       	 << fsize << " response variables" << endl;
    }
    // write out each point, one at a time
    dataItr->toFront();
    for (unsigned i = 0; i < dataItr->elementCount(); i++) {
      if (binary) {
	dataItr->currentElement().writeBinary(os);
      } else {
 	dataItr->currentElement().writeText(os);
      }
      dataItr->nextElement();
    }
  }
}

/// Read the data from a file and create a SurfDataIterator wrapper for them
void Surface::readData(std::istream& is, bool binary)
{
  delete dataItr;
  sd = new SurfData(is, binary);
  dataItr = new SurfDataIterator(*sd);
}

// ____________________________________________________________________________
// Testing
// ____________________________________________________________________________

void Surface::writeMatrix(const string header, double* mat, unsigned rows, 
  unsigned columns, ostream& os)
{
  if (header != "none" && header != "") {
    cout << header << endl;
  }
  for (unsigned r = 0; r < rows; r++) {
    for (unsigned c = 0; c < columns; c++) {
      os << setw(15) << mat[r + c*rows];
    }
    os << endl;
  }
}

void Surface::writeMatrix(const string filename, double* mat, unsigned rows, 
  unsigned columns)
{
  ofstream outfile(filename.c_str(),ios::out);
  if (!outfile) {
    cerr << "Could not open file (" << filename << ") in writeMatrix." << endl;
    return;
  }
  writeMatrix("none",mat,rows,columns,outfile);
  outfile.close();
}

#ifdef __TESTING_MODE__
  int Surface::constructCount = 0;
  int Surface::destructCount = 0;
#endif

