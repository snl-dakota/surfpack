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
#include <sstream>
#include <iomanip> 
#include "SurfPoint.h"
#include "surfpack.h"

using namespace std;
using namespace surfpack;

// ____________________________________________________________________________
// Creation, Destruction, Initialization 
// ____________________________________________________________________________

/// Initialize without any response values
SurfPoint::SurfPoint(const vector<double>& x) : x(x)
{
#ifdef __TESTING_MODE__
  constructCount++;
#endif
  init();
}

/// Initialize point with one response value
SurfPoint::SurfPoint(const vector<double>& x, double f0) : x(x), f(1)
{
#ifdef __TESTING_MODE__
  constructCount++;
#endif
  f[0] = f0;
  init();
}

/// Initialize with zero or more response values
SurfPoint::SurfPoint(const vector<double>& x, const vector<double>& f)
  : x(x), f(f)
{
#ifdef __TESTING_MODE__
  constructCount++;
#endif
  init();
}

/// Read point from istream in either text or binary format
SurfPoint::SurfPoint(unsigned xsize, unsigned fsize, istream& is, bool binary) 
{
#ifdef __TESTING_MODE__
  constructCount++;
#endif
  x.resize(xsize);
  f.resize(fsize);
  if (binary) {
    readBinary(is);
  } else { 
    readText(is);
  }
  init();
}

/// Copy constructor performs a deep copy
SurfPoint::SurfPoint(const SurfPoint& sp) : x(sp.x), f(sp.f)
{
#ifdef __TESTING_MODE__
  constructCount++;
  copyCount++;
#endif
  init();
}

/// Initialization used by all regular constructors
void SurfPoint::init()
{
  if (x.empty()) {
    throw null_point();
  }
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
  string header(
    "Error in query SurfPoint::F. Invalid responseIndex."
  );
  checkRange(header, responseIndex);
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
  string header(
    "Error in command SurfPoint::F. Invalid responseIndex. No update made."
  );
  checkRange(header, responseIndex);
  f[responseIndex] = responseValue; 
}

// ____________________________________________________________________________
// I/O
// ____________________________________________________________________________

/// Write point to stream in text or binary format
void SurfPoint::writeBinary(ostream& os) const
{
  for (unsigned i = 0; i < x.size(); i++) {
    os.write(reinterpret_cast<const char*>(&x[i]),sizeof(x[i])) ;
  }
  for (unsigned i = 0; i < f.size(); i++) {
    os.write(reinterpret_cast<const char*>(&f[i]),sizeof(f[i]));
  }
}

void SurfPoint::writeText(ostream& os) const
{
  std::_Ios_Fmtflags old_flags = os.flags();
  unsigned old_precision = os.precision(output_precision);
  os.setf(ios::scientific);
  for (unsigned i = 0; i < x.size(); i++) {
    os  << setw(field_width) << x[i] ;
  }
  for (unsigned i = 0; i < f.size(); i++) {
    os << setw(field_width) << f[i];
  }
  // Tack on some extra space before the endl to circumvent istringstream bug
  os << " " << endl;
  os.flags(old_flags);
  os.precision(old_precision);
}
 
void SurfPoint::readBinary(istream& is)
{
  try {
    // read the point in binary format
    unsigned i;
    for (i = 0; i < x.size(); i++) {
       is.read(reinterpret_cast<char*>(&x[i]),sizeof(x[i]));
    }
    for (i = 0; i < f.size(); i++) {
       is.read(reinterpret_cast<char*>(&f[i]),sizeof(f[i]));
    }
  } catch (...) {
    cerr << "Unknown error in SurfPoint::readBinary(istream& is)" << endl;
    throw;
  }
}

void SurfPoint::readText(istream& is)
{
  try {
    // read the point as text
    unsigned i;
    string sline;
    getline(is,sline);
    istringstream streamline(sline);
    for (i = 0; i < x.size(); i++) {
       streamline >> x[i];
    }
    for (i = 0; i < f.size(); i++) {
       streamline >> f[i];
    }
  } catch (...) {
    cerr << "Unknown error in SurfPoint::readText(istream& is)" << endl;
    throw;
  }
}
/// Write point to an output stream in text format
ostream& operator<<(ostream& os, const SurfPoint& sp) 
{
  sp.writeText(os);
  return os;
}

// ____________________________________________________________________________
// Testing 
// ____________________________________________________________________________

void SurfPoint::checkRange(const string& header, unsigned index) const
{
  if (index >= f.size()) {
    ostringstream errormsg;
    errormsg << header << endl;
    if (f.empty()) {
      errormsg << "There are no response values associated with this point"
               << endl;
    } else {
      errormsg << "Requested: " 
	     << index 
	     << "; actual max index: "
	     << f.size() - 1
	     << endl;
    }
    throw range_error(errormsg.str());
  }
}
#ifdef __TESTING_MODE__
  int SurfPoint::constructCount = 0;
  int SurfPoint::copyCount = 0;
  int SurfPoint::destructCount = 0;
#endif

