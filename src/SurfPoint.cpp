// Project: SURFPACK++
//
// File:       SurfPoint.cpp
// Author:     Eric Cyr 
// Modified:   Mark Richards
// 
// Description
// + SurfPoint class - a container class for a point and its response values
// + Left shift (<<) operator for SurfPoint. Easy printing.
// ____________________________________________________________________________

#include <vector>
#include <iostream>
#include <iomanip> 
#include "SurfPoint.h"

using namespace std;

// ____________________________________________________________________________
// Creation, Destruction, Initialization 
// ____________________________________________________________________________

/// Initialize without any response values
SurfPoint::SurfPoint(const vector<double>& x) : x(x)
{
#ifdef __TESTING_MODE__
  constructCount++;
#endif
}

/// Initialize point with one response value
SurfPoint::SurfPoint(const vector<double>& x, double f0) : x(x), f(1)
{
#ifdef __TESTING_MODE__
  constructCount++;
#endif
  f[0] = f0;
}

/// Initialize with zero or more response values
SurfPoint::SurfPoint(const vector<double>& x, const vector<double>& f)
  : x(x), f(f)
{
#ifdef __TESTING_MODE__
  constructCount++;
#endif
}

/// Read point from istream in either text or binary format
SurfPoint::SurfPoint(unsigned xsize, unsigned fsize, istream& is, bool binary) 
{
#ifdef __TESTING_MODE__
  constructCount++;
#endif
  unsigned i;
  x.resize(xsize);
  f.resize(fsize);
  if (!binary) {
    // read the point as text
    for (i = 0; i < xsize; i++) {
       is >> x[i];
    }
    for (i = 0; i < fsize; i++) {
       is >> f[i];
    }
  } else { 
    // read the point in binary format
    for (i = 0; i < xsize; i++) {
       is.read(reinterpret_cast<char*>(&x[i]),sizeof(x[i]));
    }
    for (i = 0; i < fsize; i++) {
       is.read(reinterpret_cast<char*>(&f[i]),sizeof(f[i]));
    }
  }
}

/// Copy constructor performs a deep copy
SurfPoint::SurfPoint(const SurfPoint& sp) : x(sp.x), f(sp.f)
{
#ifdef __TESTING_MODE__
  constructCount++;
  copyCount++;
#endif
}

/// STL data members x and f automatically cleaned up
SurfPoint::~SurfPoint() 
{
#ifdef __TESTING_MODE__
  destructCount++;
#endif
}

// ____________________________________________________________________________
// Overloaded operators 
// ____________________________________________________________________________

/// Assign sp to this unless they are already identical
SurfPoint& SurfPoint::operator=(const SurfPoint& sp) 
{
  if (*this != sp) {
    x = sp.x;
    f = sp.f;
  }
  return (*this);
}

/// Tests for deep equality
bool SurfPoint::operator==(const SurfPoint& sp) const
{
  return x == sp.x && f == sp.f;
}

/// Tests for deep inequality
bool SurfPoint::operator!=(const SurfPoint& sp) const
{
  return !(*this == sp);
} 

// ____________________________________________________________________________
// Queries 
// ____________________________________________________________________________

/// Return dimensionality of data point
unsigned SurfPoint::xSize() const
{ 
  return x.size(); 
}

/// Return number of response variables
unsigned SurfPoint::fSize() const
{ 
  return f.size(); 
}

/// Return point in the domain as an STL vector
const vector<double>& SurfPoint::X() const
{ 
  return x; 
}

/// Return response value at responseIndex
double SurfPoint::F(unsigned responseIndex) const
{ 
  if (responseIndex >= f.size()) {
    cerr << "Invalid response index. Requested: " 
	 << responseIndex 
	 << "; max: "
	 << f.size() - 1
	 << endl;
  } 
  return f[responseIndex]; 
}

// ____________________________________________________________________________
// Commands 
// ____________________________________________________________________________

/// Append a new response variable
unsigned SurfPoint::addResponse(double val)
{
  f.push_back(val);
  return f.size()-1;
}

/// Set an existing response variable to a new value
void SurfPoint::F(unsigned responseIndex, double responseValue)
{ 
  if (responseIndex >= f.size()) {
    cerr << "Invalid response index. Requested: " 
	 << responseIndex 
	 << "; max: "
	 << f.size() - 1
	 << ".  No update was made." 
	 << endl;
  } else {
    f[responseIndex] = responseValue; 
  }
}

// ____________________________________________________________________________
// I/O
// ____________________________________________________________________________

/// Write point to stream in text or binary format
ostream& SurfPoint::write(ostream& os, bool binary) const
{
  if (!binary) {
    const unsigned width = 26;
    const unsigned output_precision = 17;
    unsigned old_precision = os.precision(output_precision);
    os.setf(ios::scientific);
    for (unsigned i = 0; i < x.size(); i++) {
      os  << setw(width) << x[i] ;
    }
    for (unsigned i = 0; i < f.size(); i++) {
      os << setw(width) << f[i];
    }
    os.unsetf(ios::scientific);
    os.precision(old_precision);
  } else {
    for (unsigned i = 0; i < x.size(); i++) {
      os.write(reinterpret_cast<const char*>(&x[i]),sizeof(x[i])) ;
    }
    for (unsigned i = 0; i < f.size(); i++) {
      os.write(reinterpret_cast<const char*>(&f[i]),sizeof(f[i]));
    }
  }
  return os;
}

/// Write point to an output stream in text format
ostream& operator<<(ostream& os, const SurfPoint& sp) 
{
  return sp.write(os);
}

// ____________________________________________________________________________
// Testing 
// ____________________________________________________________________________

#ifdef __TESTING_MODE__
  int SurfPoint::constructCount = 0;
  int SurfPoint::copyCount = 0;
  int SurfPoint::destructCount = 0;
#endif

