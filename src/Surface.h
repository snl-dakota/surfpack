// These lines were automatically prepended
#ifndef SURFPACK_CONFIG_H
#define SURFPACK_CONFIG_H
#include "config.h"
#endif
// End of prepended lines

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

#include "SurfpackParser.h"

//class AbstractSurfDataIterator;
struct ErrorStruct;

/** Defines the interface for surfaces which approximate
 *  some function based on a limited number of data points.
 */
class Surface
{
protected:
class null_dimension_surface: public std::runtime_error
{
public:
  null_dimension_surface(const std::string& msg = "") 
    : std::runtime_error(msg) {}
};

class bad_metric: public std::runtime_error
{
public:
  bad_metric(const std::string& msg = "") 
    : std::runtime_error(msg) {}
};
//class too_many_terms: public std::runtime_error
//{
//public:
//  too_many_terms(const std::string& msg = "") 
//    : std::runtime_error(msg) {}
//};
  
// ____________________________________________________________________________
// Creation, Destruction, Initialization 
// ____________________________________________________________________________

protected: 
  /// Initialize SurfData to null
  Surface(); 

public:
  /// Data to be used to create surface is specified 
  Surface(SurfData* sd);

  // Makes deep copy
  Surface(const Surface& s);
  
  /// Delete data iterator 
  virtual ~Surface(); 

  /// Initialize state variables
  void init();

  /// Create a surface of the same type as 'this.'  This object's data should
  /// be replaced with the passed in, but all other attributes should
  /// be the same (e.g., a second-order polynomial should return another 
  /// second-order polynomial.  Surfaces returned by this method can be used
  /// to compute the PRESS statistic.
  virtual Surface* makeSimilarWithNewData(SurfData* sd)=0;

// ____________________________________________________________________________
// Queries 
// ____________________________________________________________________________

public:
  /// Return the name of this surface type as a string
  virtual const std::string surfaceName() const = 0;
  
  /// Return dimensionality of the surface or zero if not built
  virtual unsigned xSize();

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
  /// Make sure the surface is valid and then call evaluate
  virtual double getValue(const std::vector<double>& x); 

  /// Evaulate the approximation surface at point x and return the value.
  /// The point x must have the same dimensionality as this surface's SurfData.
  virtual double getValue(const SurfPoint& sp); 

  /// Evaluate the empirical model at the points in surfData and store 
  /// the observed and estimated values in results
//  virtual void getValue(SurfData& surfData, std::vector<double> results);

  virtual void getValue(SurfData& surfData);

  /// Evaluate the empirical model at the points in surfData and output
  /// the points and their evaluations to os.
  //virtual double test(SurfData& surfData, std::ostream& os = std::cout);
  
  /// Return the value of some error metric.
  virtual double goodnessOfFit(const std::string metricName, SurfData* surfData);
  virtual double press(SurfData& dataSet);
  virtual double rSquared(SurfData& dataSet);
  virtual double sse(SurfData& dataSet);
  virtual double mse(SurfData& dataSet);
  virtual double mrae(SurfData& dataSet);
  

// ____________________________________________________________________________
// Commands
// ____________________________________________________________________________

  /// Associates a data set with this surface object.  If this surface has
  /// already been built, it is invalidated
  virtual void setData(SurfData* sd);

  /// Set the state of the SurfData object to use the default index and points
  /// associated with this surface
  virtual void prepareData();

  /// Checks to make sure the data passed in is not null.  If it is, sets it 
  /// to point to the SurfData used to create the object.  If that is also
  /// non-existent, it is an error.  
  virtual SurfData& checkData(SurfData* dataSet);

  /// Called when the SurfData object gets modified.
  //virtual void notify();

  /// Select which response variable will be used in this surface.
  //virtual void responseIndex(unsigned index);
 
  /// Check to make sure that data are acceptable and then build.
  /// Do not build if the surface has already been built and the data have not
  /// changed
  virtual void createModel(SurfData* surfData = 0);

  /// Create an empirical model using the data from dataItr. 
  virtual void build(SurfData& data) = 0; 

  virtual void getValue(SurfData& sd, std::vector<ErrorStruct>& pts);

  virtual void config(const SurfpackParser::ArgList& arglist);

// ____________________________________________________________________________
// Helper methods 
// ____________________________________________________________________________

bool testFileExtension(const std::string& filename) const;

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
  //virtual void writeData(std::ostream& os, bool binary = false);

  /// Read the data from a file and create a SurfDataIterator wrapper for them
  //virtual void readData(std::istream& is, bool binary = false);

// ____________________________________________________________________________
// Data members 
// ____________________________________________________________________________

protected: 

  /// Number of dimensions in the data
  unsigned xsize;

  /// True if some surface has been successfully built.  This does not imply
  /// that the surface matches the current data, only that there is something
  /// valid to evaluate
  bool builtOK;

  /// True if one or more data points have been added since the surface was 
  /// built
  bool dataAdded;

  /// True if one or more data points have been modified since the surface was 
  /// built
  bool dataModified;

  /// Indices of points present in sd that were not used to make the surface 
  std::set<unsigned> excludedPoints;

  /// Data used to create this surface
  SurfData* sd;

  /// Index of the response in sd that was used to create this surface
  unsigned responseIndex;

// ____________________________________________________________________________
// Testing
// ____________________________________________________________________________

#ifdef __TESTING_MODE__
  friend class SurfaceFactoryUnitTest;
  friend class SurfaceTest;
  friend class PolynomialSurfaceTest;
 
  static int constructCount;
  static int destructCount;
#endif
};

std::ostream& operator<<(std::ostream& os,Surface&);
#endif
