// ----------------------------------------------------------
// Project: SURFPACK++
//
// File:        PolynomialSurface.h
// Author:      Mark Richards	 
// Created:     August 13, 2003
// Modified:    
//
// Description: 
// + The PolynomialSurface class does a kth order polynomial fit on the n-dimensional points
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
private:
   double average;
   PolynomialSurface() {}

public:

// constructor/destructor
////////////////////////////////

   //PolynomialSurface(const std::vector<double> &);
   PolynomialSurface(std::istream& is, int responseIndex = 0);
   PolynomialSurface(SurfData *, int order, int responseIndex = 0);
   PolynomialSurface(int numvars, int order, std::vector<double> coefficients);
   PolynomialSurface(string filename);
   //PolynomialSurface(const PolynomialSurface &);
   
   // a do nothing destructor
   //
   virtual ~PolynomialSurface(); 

// member functions
////////////////////////////////

   // get the error for this surface using the method specified
   // code: specifies the method used to calculate the error
   //virtual double getError(int code) const; 

   virtual int getMinPointCount(int dim) const;

   virtual int getDimension() const;
   //virtual ostream & writeClean(std::ostream & os = std::cout) throw(SurfException);
   virtual ostream& write(std::ostream& os = std::cout); 
   static int getMinPointCount(int order, int numvars);
   std::string getType() const;
   virtual void save(std::string filename);
    /// Save the surface to a file in binary format
   virtual void saveBinary(std::string filename);

   /// Load the surface from a file
   virtual void loadBinary(std::string filename);
   
   virtual double rSquared(AbstractSurfDataIterator* iter = 0);
   
   virtual double errorMetric(string metricName);
   virtual double press();
   virtual double rsquared();

protected:

    int k;
    int n;
    int numCoeff;
    std::vector<int> digits;
    std::vector<double> coefficients;
    int current;
   
    virtual double calculate(const std::vector<double> & x) 
      throw(SurfException);

   virtual void calculateInternals(AbstractSurfDataIterator* iter);
   void resetTermCounter();
   void nextTerm();
   double computeTerm(const std::vector<double>& x);
   void printTermLabel(std::ostream& os = std::cout);
   static int fact(int x);
   static int nChooseR(int n, int r); 
   ostream& writeMatrix(double* mat, int rows, int columns, ostream& os);
};

#endif
