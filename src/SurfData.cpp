// ----------------------------------------------------------
// Project: SURFPACK++
//
// File:        SurfData.cpp
// Author:      Eric Cyr
// Modified:	Mark Richards
//
// Description: 
// + SurfData class - this is the container for all the data 
//   from which an empirical model is created. 
// + Left shift (<<) operator for SurfData. 
// ----------------------------------------------------------

#include<vector>
#include<iostream>
#include<fstream>

#include "SurfData.h"
#include "SurfPoint.h"
//#include "Surface.h"

using namespace std;

// ____________________________________________________________________________
// Creation, Destruction, Initialization 
// ____________________________________________________________________________

/// Vector of points will be copied
SurfData::SurfData(const vector<SurfPoint>& points) 
{
#ifdef __TESTING_MODE__
  constructCount++;
#endif
  if (points.size() == 0) {
    cerr << "Error in SurfData constructor.  There are no SurfPoints" << endl;
  }
  this->xsize = points[0].xSize();
  this->fsize = points[0].fSize();
  this->points = points;
}

/// Read a set of SurfPoints from a file
SurfData::SurfData(const string filename) 
{
#ifdef __TESTING_MODE__
  constructCount++;
#endif
  // open file 
  // if filename has ".txt", read in text mode; otherwise, read in binary mode
  bool binary = (filename.find(".txt") != filename.size() - 4);
  ifstream infile(filename.c_str(), (binary ? ios::in|ios::binary : ios::in));
  if (!infile) {
    cerr << "File named \"" 
         << filename
	 << "\" could not be opened." 
	 << endl;
  }
  
  //The first data in the file should be:
  // 1. Number of points
  // 2. Dimensionality of points
  // 3. Number of response variables per point 
  unsigned size;
  if (!binary) {
    infile >> size;
    infile >> xsize;
    infile >> fsize;
  } else {
    infile.read(reinterpret_cast<char*>(&size),sizeof(size));
    infile.read(reinterpret_cast<char*>(&xsize),sizeof(xsize));
    infile.read(reinterpret_cast<char*>(&fsize),sizeof(fsize));
  }

  // Insert code for sanity check on size, xsize, and fsize

  vector< double > x(xsize);
  vector< double > f(fsize);
  unsigned index = 0;
  // loop counter
  unsigned i; 
  while (index < size && !infile.eof()); {
    if (!binary) {
      for (i = 0; i < xsize; i++) {
        infile >> x[i];
      }
      for (i = 0; i < fsize; i++) {
        infile >> f[i];
      }
    } else {
      for (i = 0; i < xsize; i++) {
        infile.read(reinterpret_cast<char*>(&x[i]),sizeof(x[i]));
      }
      for (i = 0; i < fsize; i++) {
        infile.read(reinterpret_cast<char*>(&f[i]),sizeof(f[i]));
      }
    }
    // Insert code for error checking
    SurfPoint sp(x,f);
    points.push_back(sp);
    i++;
  }
  // check to make sure end-of-file happens right here
 
  infile.close();
      
}

/// Makes a deep copy of the object 
SurfData::SurfData(const SurfData& sd) 
{
#ifdef __TESTING_MODE__
  constructCount++;
  copyCount++;
#endif
  if (sd.size() == 0) {
    cerr << "Unacceptable copy.  SurfData to be copied has no data." << endl;
  }
  this->xsize = sd.xsize;
  this->fsize = sd.fsize;
  this->points = sd.points;
}

/// STL data members' resources automatically deallocated 
SurfData::~SurfData() 
{
  // clean up the related surfaces
  //listeners.clear();
  // clean up all the created points
  //deletePoints();
#ifdef __TESTING_MODE__
  destructCount++;
#endif
}

// ____________________________________________________________________________
// Overloaded operators 
// ____________________________________________________________________________

/// makes deep copy 
SurfData& SurfData::operator=(const SurfData& sd)
{
  if (*this != sd) {
    this->xsize = sd.xsize;
    this->fsize = sd.fsize;
    this->points = sd.points;
    // now all the surfaces need to be update
    // or at least told of a change
    //notifyListeners();
  }
  return (*this);
}

/// makes deep comparison
bool SurfData::operator==(const SurfData& sd)
{
  return (this->xsize == sd.xsize && 
          this->fsize == sd.fsize &&
          this->points == sd.points);
}
      
/// makes deep comparison
bool SurfData::operator!=(const SurfData& sd)
{
  return !(*this == sd);
}

// ____________________________________________________________________________
// Queries 
// ____________________________________________________________________________

/// Return the number of SurfPoints in the data set 
unsigned SurfData::size() const 
{ 
  return points.size(); 
}

/// Return the dimensionality of the SurfPoints 
unsigned SurfData::xSize() const 
{ 
  return xsize; 
}

/// Return the number of response values in the data set
unsigned SurfData::fSize() const 
{ 
  return fsize; 
}

/// Return a point from the data set
SurfPoint& SurfData::Point(unsigned index) 
{
  // What if there are no points?
  if (points.size() == 0) {
    cerr << "Cannot return points[" << index << "].  There are no points." << endl;
  }
  if (index >= points.size()) {
    cerr << "Out of range in SurfData::Point; returning last point" << endl;
    return points[points.size()-1];
  }
  return points[index];
}

/// Return a reference to the SurfPoints vector 
vector<SurfPoint>& SurfData::Points() 
{ 
  return points; 
}

// ____________________________________________________________________________
// Commands 
// ____________________________________________________________________________

/// Add a point to the data set. The parameter point will be copied.
void SurfData::addPoint(const SurfPoint& sp) 
{
  points.push_back(sp);
}

/// Add a new response variable to each point. 
/// Return the index of the new variable.
unsigned SurfData::addResponse()
{
  unsigned newindex;
  if (points.size() == 0) {
    cerr << "Cannot add another response because there are no points." << endl;
    return 0;
  } else {
    newindex = points[0].addResponse(0);
    fsize = newindex + 1;
    for (unsigned i = 1; i < points.size(); i++) {
      if (points[i].addResponse() != newindex) {
        cerr << "Size mismatch among data points in this SurfData object."
             << endl
             << "This will likely cause a fatal error." 
             << endl;
      }
    }
  }
  return newindex;
}

//void SurfData::addListener(Surface * surface)
//{
//  /// only add the listener if its not already there
//  list<Surface*>::iterator itr = 
//    find(listeners.begin(),listeners.end(),surface);
//  if(itr==listeners.end()) {
//    listeners.push_back(surface);
//  }
//}
//
//void SurfData::removeListener(Surface * surface)
//{
//  // make sure its OK to erase the object, then do so
//  list<Surface*>::iterator itr =
//    find(listeners.begin(),listeners.end(),surface);
//  if(itr!=listeners.end()) {
//    listeners.erase(itr);
//  }
//}

// ____________________________________________________________________________
// I/O
// ____________________________________________________________________________

/// Write a set of SurfPoints to a file.  Opens the file and calls other version
void SurfData::write(const std::string filename) const
{
  bool binary = (filename.find(".txt") != filename.size() - 4);
  ofstream outfile(filename.c_str(), (binary ? ios::out|ios::binary : ios::out));
  if (!outfile) {
    cerr << "Could not open " << filename << " for writing." << endl;
  } else {
    write(outfile, binary);
  }
}

/// Write a set of SurfPoints to an output stream
ostream& SurfData::write(ostream& os, bool binary) const
{
  if (!binary) {
    os << points.size() << endl
       << xsize << endl 
       << fsize << endl;
    os << 0 << endl; // grad size
  } else {
    unsigned s = points.size();
    os.write((char*)&s,sizeof(s));
    os.write((char*)&xsize,sizeof(xsize));
    os.write((char*)&fsize,sizeof(fsize));
  }
  vector<SurfPoint>::const_iterator itr;
  itr = points.begin();
  while (itr != points.end()) {
    (*itr).write(os,binary);
    if (!binary) {
      // Surfpoint->write(os) does not write newline after each point
      os << endl; 
    } 
    ++itr;
  }
  return os;
}

// so a SurfData object can be printed
ostream& operator<<(ostream& os, const SurfData& sd) 
{ 
  return sd.write(os); 
}

// ____________________________________________________________________________
// Helper methods 
// ____________________________________________________________________________

//void SurfData::notifyListeners()
//{
//  list<Surface*>::iterator itr;
//  for(itr=listeners.begin();itr!=listeners.end();++itr)
//    (*itr)->notify();
//}

// ____________________________________________________________________________
// Testing 
// ____________________________________________________________________________

#ifdef __TESTING_MODE__
  int SurfData::constructCount = 0;
  int SurfData::copyCount = 0;
  int SurfData::destructCount = 0;
#endif

