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
#include "Surface.h"
#include "SurfPoint.h"
#include "SurfData.h"
#include "SurfDataIterator.h"


using namespace std;

// ____________________________________________________________________________
// Creation, Destruction, Initialization 
// ____________________________________________________________________________

Surface::Surface() : valid(false), originalData(false), dataItr(0)
{
#ifdef __TESTING_MODE__
  constructCount++;
#endif
}

Surface::Surface(AbstractSurfDataIterator* dataItr) :
  valid(false), originalData(false), dataItr(dataItr) 
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
    currentPoint.write(outfile);
    outfile << setw(20) << estimatedResponse 
            << setw(20) << (currentPoint.F() - estimatedResponse)
  	  << endl;	
    error += (currentPoint.F() - estimatedResponse) * 
      (currentPoint.F() - estimatedResponse);
  }
  outfile.close();
  return error;
}    

double Surface::errorMetric(string metricName)
{
  if (metricName == "press") {
    cout << "Call Press from Surface" << endl;
    return press();
  } else {
    cout << "No error metric of that type in this class" << endl;
  }
  return 0;
}

double Surface::press()
{
  cout << "Surface::press" << endl;
  return 0;
}

double Surface::rSquared()
{
  cout << "Surface::rsquared" << endl;
  return 0;
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

// ____________________________________________________________________________
// I/O 
// ____________________________________________________________________________

void Surface::write(const string filename)
{
  bool binary = (filename.find(".txt") != filename.size() - 4);
  ofstream outfile(filename.c_str(), (binary ? ios::out|ios::binary : ios::out));
  if (!outfile) {
    cerr << "File named \"" 
         << filename
	 << "\" could not be opened." 
	 << endl;
    return;
  } else if (binary) {
    writeBinary(outfile);
  } else {
    writeText(outfile);
  }
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
    readBinary(infile);
  } else {
    readBinary(infile);
  }
  infile.close();
}

ostream& operator<<(ostream& os,Surface& surface)
{ 
  surface.writeText(os); 
  return os;
}


// ____________________________________________________________________________
// Testing
// ____________________________________________________________________________

#ifdef __TESTING_MODE__
  int Surface::constructCount = 0;
  int Surface::destructCount = 0;
#endif

