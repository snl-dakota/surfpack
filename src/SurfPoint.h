// These lines were automatically prepended
#ifndef SURFPACK_CONFIG_H
#define SURFPACK_CONFIG_H
#include "config.h"
#endif
// End of prepended lines

// Project: SURFPACK++
//
// File:       SurfPoint.h
// Author:     Eric Cyr 
// Modified:   Mark Richards
// 
// Description
// + SurfPoint class - a container class for a point and its response values
// + Left shift (<<) operator for SurfPoint. Easy printing.
// ____________________________________________________________________________

#ifndef __SURF_POINT_H__
#define __SURF_POINT_H__

// INVARIANTS: for SurfPoint
// ---------------------------

// not yet specified
#include <stdexcept>


class SurfPoint {

// Nested Exception class used when an attempt is made to create a SurfPoint 
// with 0 dimensions
private:
class null_point : public std::runtime_error
{
public:
  null_point(const std::string& msg = 
    "Error: attempt to make SurfPoint with 0 dimensions.") 
    : std::runtime_error(msg) {}
};
  
// ____________________________________________________________________________
// Creation, Destruction, Initialization 
// ____________________________________________________________________________
public:

  /// Initialize without any response values
  SurfPoint(const std::vector<double>& x);

  /// Initialize point with one response value
  SurfPoint(const std::vector<double>& x, double f0);
  
  /// Initialize with zero or more response values
  SurfPoint(const std::vector<double>& x, const std::vector<double>& f);
  
  /// Read point from istream in either text or binary format
  SurfPoint(unsigned xsize, unsigned fsize, std::istream& is, bool binary = false);

  /// Copy constructor performs a deep copy
  SurfPoint(const SurfPoint& sp);

  /// STL data members x and f automatically cleaned up
  ~SurfPoint();

private:
  /// Initialization used by all regular constructors
  void init();

protected:
  /// Default constructor explicitly disallowed.  A SurfPoint must at least
  /// specify a point in some domain space.
  SurfPoint();

// ____________________________________________________________________________
// Overloaded operators 
// ____________________________________________________________________________
public:

  /// Assign sp to this unless they are already identical
  SurfPoint& operator=(const SurfPoint& sp);

  /// Tests for deep equality
  bool operator==(const SurfPoint& sp) const;
  
  /// Tests for deep inequality
  bool operator!=(const SurfPoint& sp) const;
  
  /// Function object for use with sets of SurfPoint objects (in particular,
  /// a SurfData object)
  class SurfPointPtrLessThan
  {
  public:
    bool operator()(const SurfPoint* sp1, const SurfPoint* sp2);
  };
      

// ____________________________________________________________________________
// Queries 
// ____________________________________________________________________________

  /// Return dimensionality of data point
  unsigned xSize() const;

  /// Return number of response variables
  unsigned fSize() const;

  /// Return point in the domain as an STL vector
  const std::vector<double>& X() const;

  /// Return response value at responseIndex
  double F(unsigned responseIndex = 0) const;

// ____________________________________________________________________________
// Commands 
// ____________________________________________________________________________

  /// Append a new response variable
  unsigned addResponse(double val = 0); 

  /// Set an existing response variable to a new value
  void F(unsigned responseIndex, double responseValue);
// ____________________________________________________________________________
// I/O
// ____________________________________________________________________________

  /// Write point to stream in text or binary format
  void writeBinary(std::ostream& os) const;
  void writeText(std::ostream& os) const;
  void readBinary(std::istream& is);
  void readText(std::istream& is);

// ____________________________________________________________________________
// Data members 
// ____________________________________________________________________________

protected:
  /// The point in the domain.  The size of the vector
  /// is the dimensionality of the space. 
  std::vector<double> x;          

  /// Zero or more response values at x (i.e., f1(x), f2(x) ... )
  std::vector<double> f;      

// ____________________________________________________________________________
// Testing 
// ____________________________________________________________________________
protected:

void checkRange(const std::string& header, unsigned index) const;


#ifdef __TESTING_MODE__
  friend class SurfPointTest;
  friend class SurfDataTest;
public:
  static int constructCount;
  static int copyCount;
  static int destructCount;
#endif

};

/// Write point to an output stream in text format
std::ostream& operator<<(std::ostream& os, const SurfPoint& sp); 

#endif
