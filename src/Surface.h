// ----------------------------------------------------------
// Project: SURFPACK++
//
// File:        Surface.h
// Author:      Tony Giunta
// Created:     June 13, 2002
// Modified:    Eric Cyr - June 18, 2002        
//                         June 20, 2002
//                         June 26, 2002
//                         July 15, 2002
// Modified:	Mark Richards - June 23, 2003
// 				August 12, 2003
//
// Description: 
// + The Surface class provides an interface for the
//   different surface types, this class cannot be instantiated
// ----------------------------------------------------------

#ifndef __SURFACE_H__
#define __SURFACE_H__

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <cassert>

#include "SurfException.h"
#include "SurfPoint.h"
#include "SurfDataIterator.h"
#include "SkipSurfDataIterator.h"



class SurfData;

/** Defines the interface for surfaces which approximate
 *  some function based on a limited number of data points.
 */
class Surface
{
public:
static const int polynomialSurfaceID = 1;
static const int krigingSurfaceID = 2;

// Constructor/destructor
////////////////////////////////
   
   /// A do-nothing destructor
   virtual ~Surface(); 

// member functions
////////////////////////////////

   /// Select which response variable will be used in this surface
   virtual bool setResponseIndex(int);

   /// Return the index of the response variable used in this surface
   virtual int getResponseIndex();
 
   /// Return the value of the PRESS error metric
   //virtual double getPressStatistic();

   /// Return the value of some error metric
   virtual double errorMetric(string metricName);
   virtual double press();
   virtual double rsquared();

   /// Return the minumum number of points needed to create
   ///     a surface of this type 
   virtual int getMinPointCount() const
   { return getMinPointCount(getDimension()); }

   /// Return the minimum number of points needed to create
   /// a surface of this type for a function with dim dimensions
   virtual int getMinPointCount(int dim) const = 0;
 
   /// Return the number of dimensions in the SurfData used to make this surface
   virtual int getDimension() const = 0;

   /// Set the iterator for this class. 
   /// Uses the clone funtion to copy the object.
   /// Copy constructors are not polymorphic
   //void setIterator(const AbstractDataIterator & iter);
  
   /// Set the points from which this surface will be calculated
   //virtual void setData(SurfData * data);

   /// Return the SurfData object this surface is using
   virtual SurfData * getData() const { return surfData; }

   /// Save the surface to a file
   virtual void save(std::string filename) = 0;

   /// Save the surface to a file in binary format
   virtual void saveBinary(std::string filename);

   /// Load the surface from a file
   virtual void loadBinary(std::string filename);
   
   
   /// Return true if there is a SurfData object associated with this class
   //virtual bool isData() { return surfData!=0; }

   /// Return the "error" for this surface using the metric specified by code
   //virtual double getError(int code) const = 0; 

   /// Specify the coefficients that encode this surface (as opposed to 
   /// computing them with Surface::calculateInternals() 
   //virtual void setCoefficients(const std::vector<double> & coeff) 
   //{ coefficients = coeff; needsRebuild = true;  }

   /// Retrieve the coefficients that encode this surface
   //virtual void getCoefficients(std::vector<double> & coeff)
   //{ if(needsRebuild) build(); coeff = coefficients; }

   /// Evaulate the approximation surface at point x and return the value.
   /// The point x must have the same dimensionality as this surface's SurfData object
   virtual double evaluate(const std::vector<double> & x) throw(SurfException);

   /// Evaluate the empirical model at the points in surfData and output
   /// the points and their evaluations to os
   virtual void evaluate(SurfData* surfData);

   /// Evaluate the empirical model at the points in surfData and output
   /// the points and their evaluations to os
   virtual double test(SurfData* surfData, std::ostream & os = std::cout);
   
   // evaulate the surface at a given point and return the value
   // simply calls the other evaluate (only one to overload!!!)
   //
   // x: the point to evaluate the surface
   // dim: the # of dimensions in the point 
   // return: returns the value at the point x
   //
   //virtual double evaluate(double * x,int dim) throw(SurfException);

   // evaulate the surface at a given point and return the value
   // and gradient
   //
   // x: the point to evaluate the surface
   // grad: where the gradient is returned
   // return: returns the value at the point x
   //
   //virtual double evaluate(const std::vector<double> & x,std::vector<double> & grad)
   //   throw(SurfException);

   // evaulate the surface at a given point and return the value
   // and gradient
   //
   // x: the point to evaluate the surface
   // grad: where the gradient is returned
   // dim: the # of dimensions in the point 
   // return: returns the value at the point x
   //
   //virtual double evaluate(const double * x,double * grad,int dim) throw(SurfException);

   // evaluate the gradient of the surface at a given point
   //
   // x: the point ot evaluate the gradient at
   // grad: the location where the gradient will be
   //
   //virtual void evaluateGrad(const std::vector<double> & x,std::vector<double> & grad)
   //   throw(SurfException);

   // evaluate the gradient of the surface at a given point
   //
   // x: the point ot evaluate the gradient at
   // grad: the location where the gradient will be
   // xGradSz: the size of the x and the grad arrays
   //
   //virtual void evaluateGrad(const double * x,double * grad,int xGradSz)
   //   throw(SurfException);

   /// Rebuild the surface using the points in the SurfData object 
   virtual void build() throw(SurfException);
   virtual void build(AbstractSurfDataIterator* iter) throw(SurfException);

   /// Return the name of this surface type as a string
   virtual std::string getType() const = 0;
   //{ return surfType; }
   
   /// Print the surface's information to the output stream os.
   /// The stream is returned to enable cascading. 
   virtual std::ostream& write(std::ostream & os = cout) = 0;
//   virtual std::istream& read(std::istream & is = cin) = 0;

   /// Print the surface's information to the output stream os.
   /// The stream is returned to enable cascading. 
   /// The output is neater and more concise than with print.
   //virtual std::ostream & printClean(std::ostream & os = cout) throw(SurfException);

   /// Called when the SurfData object gets modified
   virtual void notify()
   { if(surfData!=0) needsRebuild = true; }

   /// Compute and return the value of the PRESS error metric.
   virtual double computePressStatistic();

#ifdef __TESTING_MODE__
   friend class SurfaceFactoryUnitTest;
   friend class SurfaceUnitTest;
  
   static int constructCount;
   static int destructCount;
#endif

protected:

   // These constructors are here (protected) to "hide"
   // them from the user
   //Surface(const std::string &);
   //Surface(const std::vector<double> &);
   Surface(SurfData * surfData, int responseCount = 0);
   //Surface(const std::string &,SurfData *,const std::vector<double> &);
   Surface(const Surface & s);

   // Evaluate the approximation surface at point x, (including gradient). 
   // It is called by evaluate
   //virtual double calculate(const std::vector<double> & x,
   //                         std::vector<double> & grad,
   //                         bool isGrad=true) 
   //   throw(SurfException) = 0;

   // Evalute the approximation surface at point x.
   // Gradient is not computed.
   // Called by evaluate(x).
   virtual double calculate(const std::vector<double> & x) 
      throw(SurfException) = 0;

   // Use the SurfData object to create an approximation surface
   // Each child class implements its own method of approximation.
   virtual void calculateInternals(AbstractSurfDataIterator* iter) = 0;

//   std::string surfType;              // the name of this surface type
 
   SurfData * surfData;               // the data object to use
                                      // when building the surface

 //  std::vector<double> coefficients;  // the coefficients that represent
                                      // this surface

   bool needsRebuild;                 // do we need to recaculate this surface

   int responseIndex;                 // which surface do we use

  // AbstractDataIterator * surfDataIter;   // iterator to used for SurfData

   //double pressStatistic;		// "leave one out" error metric for surface fit


   // the default constructor is not used
   Surface() {}
};

std::ostream & operator<<(std::ostream & os,Surface &);

// allows for a tolerance to be set which helps determines how
// small a coefficient has to be before it simply prints out as zero
// in the printClean method
bool significant(double val);
#endif
