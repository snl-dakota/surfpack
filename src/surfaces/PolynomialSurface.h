// These lines were automatically prepended
#ifndef SURFPACK_CONFIG_H
#define SURFPACK_CONFIG_H
#include "config.h"
#endif
// End of prepended lines

// ----------------------------------------------------------
// Project: SURFPACK++
//
// File:        PolynomialSurface.h
// Author:      Mark Richards	 
// Created:     August 13, 2003
// Modified:    
//
// Description: 
// + The PolynomialSurface class does a kth order polynomial fit 
//   on the n-dimensional points
// ----------------------------------------------------------

#ifndef __POLYNOMIAL_SURFACE_H__
#define __POLYNOMIAL_SURFACE_H__

#include <iostream>
#include <vector>
#include <string>
#include <stdexcept>
#include "SurfpackParser.h"

class SurfPoint;
class Surface;

class PolynomialSurface : public Surface
{

// ____________________________________________________________________________
// Creation, Destruction, Initialization 
// ____________________________________________________________________________

public:
  
  PolynomialSurface(SurfData* sd, unsigned order = 2);
  PolynomialSurface(unsigned xsize, unsigned order, 
    std::vector<double> coefficients);
  PolynomialSurface(const std::string filename);
  PolynomialSurface(const PolynomialSurface& other);
  virtual void config(const SurfpackParser::Arg& arg);
  
  /// Create a surface of the same type as 'this.'  This objects data should
  /// be replaced with the dataItr passed in, but all other attributes should
  /// be the same (e.g., a second-order polynomial should return another 
  /// second-order polynomial.  Surfaces returned by this method can be used
  /// to compute the PRESS statistic.
  virtual PolynomialSurface* makeSimilarWithNewData(SurfData* surfData);
 
  virtual ~PolynomialSurface(); 

private:
   /// explicitly disallow default constructor
   PolynomialSurface(); 

// ____________________________________________________________________________
// Queries 
// ____________________________________________________________________________

public:

  virtual const std::string surfaceName() const;
  static unsigned minPointsRequired(unsigned xsize, unsigned order);
  virtual unsigned minPointsRequired() const;
  virtual double evaluate(const std::vector<double>& x); 
  //virtual double errorMetric(std::string metricName);
  //virtual double press();
  //virtual double rSquared(AbstractSurfDataIterator* iter = 0);

// ____________________________________________________________________________
// Commands 
// ____________________________________________________________________________

  virtual void build(SurfData& data);


// ____________________________________________________________________________
// Helper methods 
// ____________________________________________________________________________

  //static unsigned fact(unsigned x);
  //static unsigned nChooseR(unsigned n, unsigned r); 
  void resetTermCounter() const;
  double computeTerm(const std::vector<double>& x) const;
  void nextTerm() const;
// ____________________________________________________________________________
// Data members 
// ____________________________________________________________________________

protected:
  static const std::string name;
  unsigned order;
  std::vector<double> coefficients;
  mutable std::vector<unsigned> digits;
public:
  mutable unsigned termIndex;
  mutable bool lastTerm;
//protected:
   
// ____________________________________________________________________________
// I/O 
// ____________________________________________________________________________

   /// Save the surface to a file in binary format
  virtual void writeBinary(std::ostream& os);

  virtual void writeText(std::ostream& os);

  /// Load the surface from a file
  virtual void readBinary(std::istream& is);

  virtual void readText(std::istream& is);

  virtual void printTermLabel(std::ostream& os = std::cout);
  virtual void printTermComponents(std::ostream& os = std::cout);

// ____________________________________________________________________________
// Testing 
// ____________________________________________________________________________

#ifdef __TESTING_MODE__ 
  friend class PolynomialSurfaceTest;

public:
  static int constructCount;
  static int copyCount;
  static int destructCount;
#endif
};

#endif
