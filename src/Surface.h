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
struct ErrorStruct;

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
  virtual double evaluate(const std::vector<double>& x) = 0; 

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
  virtual double errorMetric(const std::string metricName, 
    AbstractSurfDataIterator* itr);
  virtual double press(AbstractSurfDataIterator* itr);
  virtual double rSquared(AbstractSurfDataIterator* itr);
  virtual double sse(AbstractSurfDataIterator* itr);
  virtual double mse(AbstractSurfDataIterator* itr);
  

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

  /// Create a surface of the same type as 'this.'  This objects data should
  /// be replaced with the dataItr passed in, but all other attributes should
  /// be the same (e.g., a second-order polynomial should return another 
  /// second-order polynomial.  Surfaces returned by this method can be used
  /// to compute the PRESS statistic.
  virtual Surface* makeSimilarWithNewData(AbstractSurfDataIterator* dataItr)=0;

  virtual void evaluate(AbstractSurfDataIterator* itr, 
    std::vector<ErrorStruct>& pts);

// ____________________________________________________________________________
// I/O 
// ____________________________________________________________________________

  /// Write the surface out to a file.  Files with extension .txt are written
  /// out in text mode; others are written in binary mode.
  virtual void write(const std::string filename);

  /// Read the surface from a file.  Files with extension .txt are read in text 
  /// mode; others are read in binary mode.
  virtual void read(const std::string filename);

  /// Write the surface in binary format
  virtual void writeBinary(std::ostream& os) = 0; 

  /// Write the surface in text format
  virtual void writeText(std::ostream& os) = 0; 

  /// Read the surface in binary format
  virtual void readBinary(std::istream& is) = 0; 

  /// Read the surface in text format
  virtual void readText(std::istream& is) = 0; 

  /// Write the associated data to a stream.  Not all iterators will use all of 
  /// the data available in their SurfData objects, so the writing must 
  /// necessarily go through the iterator.
  virtual void writeData(std::ostream& os, bool binary = false);

  /// Read the data from a file and create a SurfDataIterator wrapper for them
  virtual void readData(std::istream& is, bool binary = false);

// ____________________________________________________________________________
// Data members 
// ____________________________________________________________________________

protected: 

  /// Is this surface in a valid, consistent state? 
  bool valid;                 

  /// Set to true when data used to build surface is present 
  bool originalData;

  /// Collection of points used to create empirical model
  AbstractSurfDataIterator* dataItr;               

  /// Created and passed to dataItr only if surface is created from a file
  SurfData* sd;

// ____________________________________________________________________________
// Testing
// ____________________________________________________________________________

void writeMatrix(const std::string header, double* mat, unsigned rows, 
  unsigned columns, std::ostream& os);
void writeMatrix(const std::string filename, double* mat, unsigned rows, 
  unsigned columns);
#ifdef __TESTING_MODE__
  friend class SurfaceFactoryUnitTest;
  friend class SurfaceUnitTest;
 
  static int constructCount;
  static int destructCount;
#endif
};

std::ostream& operator<<(std::ostream& os,Surface&);
#endif
