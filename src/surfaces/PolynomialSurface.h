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

#include "AbstractSurfDataIterator.h"
#include "SurfException.h"
#include "SurfPoint.h"
#include "Surface.h"

class PolynomialSurface : public Surface
{

// ____________________________________________________________________________
// Creation, Destruction, Initialization 
// ____________________________________________________________________________

public:
  
  PolynomialSurface(SurfData& sd, unsigned order, unsigned responseIndex = 0);
  PolynomialSurface(AbstractSurfDataIterator* dataItr, unsigned order);
  PolynomialSurface(unsigned xsize, unsigned order, 
    std::vector<double> coefficients);
  PolynomialSurface(const std::string filename);
  //PolynomialSurface(const PolynomialSurface&);
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

  virtual void build();

  /// Create a surface of the same type as 'this.'  This objects data should
  /// be replaced with the dataItr passed in, but all other attributes should
  /// be the same (e.g., a second-order polynomial should return another 
  /// second-order polynomial.  Surfaces returned by this method can be used
  /// to compute the PRESS statistic.
  virtual PolynomialSurface* makeSimilarWithNewData
    (AbstractSurfDataIterator* dataItr);

// ____________________________________________________________________________
// Helper methods 
// ____________________________________________________________________________

  static unsigned fact(unsigned x);
  static unsigned nChooseR(unsigned n, unsigned r); 
  void resetTermCounter();
  double computeTerm(const std::vector<double>& x);
  void nextTerm();
// ____________________________________________________________________________
// Data members 
// ____________________________________________________________________________

protected:
  static const std::string name;
  unsigned xsize;
  unsigned order;
  
  std::vector<double> coefficients;
  std::vector<unsigned> digits;
  unsigned termIndex;
   
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

// ____________________________________________________________________________
// Testing 
// ____________________________________________________________________________

#ifdef __TESTING_MODE__ 
  friend class SurfDataUnitTest;
  friend class SurfaceUnitTest;

public:
  static int constructCount;
  static int copyCount;
  static int destructCount;
#endif
};

#endif
