// ----------------------------------------------------------
// Project: SURFPACK++
//
// File:        Surface.h
// Author:      Tony Giunta
// Modified:    Eric Cyr 
// Modified:	Mark Richards 
//
// Description: 
// + The Surface class provides an interface for the
//   different surface types; this class cannot be instantiated
// ----------------------------------------------------------

#ifndef __SURFACE_H__
#define __SURFACE_H__

class SurfPoint;
class SurfData;
class AbstractSurfDataIterator;

/** Defines the interface for surfaces which approximate
 *  some function based on a limited number of data points.
 */
class Surface
{
// ____________________________________________________________________________
// Creation, Destruction, Initialization 
// ____________________________________________________________________________

public:
   
  /// Initialize dataItr to null
  Surface(); 

  // Iterator wraps data that will be used to construct the surface
  Surface(AbstractSurfDataIterator* dataItr);

  // Makes deep copy
  Surface(const Surface& s);
  
  /// Delete data iterator 
  virtual ~Surface(); 

// ____________________________________________________________________________
// Queries 
// ____________________________________________________________________________

public:
  /// Return the name of this surface type as a string
  virtual const std::string surfaceName() const = 0;
  
  /// Return true if surface has been properly built and is in a consistent
  /// state
  bool isValid() const;

  /// Return true if the data that was used to create the surface is available.
  /// Some error metrics require the original data; others do not.
  bool hasOriginalData() const;

  /// Return true if there is a sufficient number of correctly formatted data
  /// points
  bool acceptableData() const;

  /// Return the minumum number of points needed to create s surface of this
  /// type. 
  virtual unsigned minPointsRequired() const = 0;
 
  /// Return the index of the response variable used in this surface.
  //virtual unsigned responseIndex();
 
  /// Evaulate the approximation surface at point x and return the value.
  /// The point x must have the same dimensionality as this surface's SurfData.
  virtual double evaluate(const std::vector<double>& x); 

  /// Evaulate the approximation surface at point x and return the value.
  /// The point x must have the same dimensionality as this surface's SurfData.
  virtual double evaluate(const SurfPoint& sp); 

  /// Evaluate the empirical model at the points in surfData and output
  /// the points and their evaluations to os.
  virtual void evaluate(SurfData& surfData);

  /// Evaluate the empirical model at the points in surfData and output
  /// the points and their evaluations to os.
  virtual double test(SurfData& surfData, std::ostream& os = std::cout);
  
  /// Return the value of some error metric.
  virtual double errorMetric(std::string metricName);
  virtual double press();
  virtual double rSquared();

// ____________________________________________________________________________
// Commands
// ____________________________________________________________________________

  /// Called when the SurfData object gets modified.
  //virtual void notify();

  /// Select which response variable will be used in this surface.
  //virtual void responseIndex(unsigned index);
 
  /// Rebuild the surface, if necessary
  virtual void ensureValidity();

  /// Create an empirical model using the data from dataItr. 
  virtual void build() = 0; 

// ____________________________________________________________________________
// I/O 
// ____________________________________________________________________________

  /// Write the surface out to a file.  Files with extension .txt are written
  /// out in text mode; others are written in binary mode.
  virtual void write(std::string filename);

  /// Read the surface from a file.  Files with extension .txt are read in text 
  /// mode; others are read in binary mode.
  virtual void read(std::string filename);

  /// Write the surface in binary format
  virtual void writeBinary(std::ostream& os) = 0; 

  /// Write the surface in text format
  virtual void writeText(std::ostream& os) = 0; 

  /// Read the surface in binary format
  virtual void readBinary(std::istream& is) = 0; 

  /// Read the surface in text format
  virtual void readText(std::istream& is) = 0; 

// ____________________________________________________________________________
// Data members 
// ____________________________________________________________________________

protected: 

  // Is this surface in a valid, consistent state? 
  bool valid;                 

  // Set to true when data used to build surface is present 
  bool originalData;

  // Collection of points used to create empirical model
  AbstractSurfDataIterator* dataItr;               

// ____________________________________________________________________________
// Testing
// ____________________________________________________________________________

#ifdef __TESTING_MODE__
  friend class SurfaceFactoryUnitTest;
  friend class SurfaceUnitTest;
 
  static int constructCount;
  static int destructCount;
#endif
};

std::ostream& operator<<(std::ostream& os,Surface&);
#endif
