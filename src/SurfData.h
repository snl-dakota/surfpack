// ----------------------------------------------------------
// Project: SURFPACK++
//
// File:        SurfData.h
// Author:      Eric Cyr
// Created:     June 6, 2002
// Modified:    June 17, 2002
//              June 24 - 26, 2002
//              July 11, 2002
// Modified:	Mark Richards
// 		June 24, 2003
// 		August 15, 2003
//
// Description: 
// + SurfData class - this is the container for all the data 
//   from which a surrogate/model is created. 
// + Left shift (<<) operator for SurfData. 
// ----------------------------------------------------------

#ifndef __SURF_DATA_H__
#define __SURF_DATA_H__

#include <iostream>
#include <vector>
#include <string>
#include <list>
#include <algorithm>

#include "SurfException.h"
#include "SurfPoint.h"

class Surface;

class SurfData
{
private:

// private data memebers
////////////////////////////////////
  
   /// Dimensionality of the space from wich the SurfPoints are drawn
   int dimension;

   /// Number of response variables in the data set 
   int responses;

   /// The set of points in this data set
   std::vector<SurfPoint*> points; 
 
   /// List of pointers to listening/observing Surface objects 
   /// which need to be notified when this object changes
   std::list<Surface*> listeners;

// private helper methods 
////////////////////////////////////
 
   /// Notify listening surfaces whenever a SurfPoint
   /// is added or removed.
   virtual void notifyListeners(); 

   /// Return the number of SurfPoints in the data set 
   //virtual int getPointCount() const;
   


public:

#ifdef __TESTING_MODE__ 
   friend class SurfDataUnitTest;
   friend class SurfaceUnitTest;
 
   static int constructCount;
   static int destructCount;
#endif

// constructors/destructor
////////////////////////////////////
 
   // Constructs a SurfData object
   SurfData();

   /// copy constructor
   SurfData(const SurfData & data); 

   // cleans up after a SurfData object
   virtual ~SurfData();

// public methods
////////////////////////////////////

   /// Return the number of SurfPoints in the data set 
   virtual int size() const;
   
   /// retrieve a point from the data set
   SurfPoint* getPoint(int x);

   /// Delete all the points in the "points" 
   /// vector
   virtual void deletePoints();

   /// Remove all the points in "points"
   /// but do not delete them (in the sense of deallocating the memory)
   //virtual void removePoints();

   /// Add a new response variable for each point. 
   /// Return the index of the new variable.
   int addResponse(); 

   /// Make a "shallow" copy of the point. That
   /// is only copy the pointer to the point
   /// and not the entire data structure
   virtual void shallowAddPoint(SurfPoint * pt);
   
   /// Add a point to the data set. The parameter point will be copied
   virtual void addPoint(const SurfPoint & point);

   // add a point to the data set
   //
   // point: a vector to the point to be added
   // fval: the function value at this point
   //
   //virtual void addPoint(const std::vector<double> & point,double fval);

//   virtual void addPoint(const std::vector<double> & point,
//                         double fval,
//                         const std::vector<double> & grad);

   // add a point to the data set
   //
   // point: a vector to the point to be added
   // fvals: the response values at this point
   //
   //virtual void addPoint(const std::vector<double> & point,
    //                     const std::vector<double> & fvals);

   // add a point to the data set
   //
   // point: a vector to the point to be added
   // fvals: the response values at this point
   // grads: the gradients at this point
   //
//   virtual void addPoint(const std::vector<double> & point,
//                         const std::vector<double> & fvals,
//                         const std::vector<std::vector<double> > & grads);

   // add a point to the data set
   // 
   // point: an array of "dim" values representing the domain
   //        of the point to be added
   // fvals: an array of "respCnt" values representing the response
   //        values at this point
   // grads: an array of "dim * respCnt" values representing the
   //        gradient values at this point
   // dim: dimension of point (the dimension of the domain)
   // respCnt: the number of responses in this point
   //
   // 
//   virtual void addPoint(const double * point,
//                         const double * fvals, 
//                         const double * grads,
//                         int dim,int respCnt);

   // add a set of points to the data
   // 
   // points: the vector of points to be added 
   // 
//   virtual void addPoints(const std::vector<SurfPoint*> & points);

   // the "FORTRAN" style interface to the SurfData object
   //
   // domain: the x values
   // response: the various response values
   //
 //  virtual void addPoints(const std::vector<std::vector<double> > & domain,
 //                         const std::vector<std::vector<double> > & response) 
 //     throw(SurfException);

   // the "FORTRAN" style interface to the SurfData object
   //
   // domain: the x values
   // response: the various response values
   // grad: the various response gradients
   //
 //  virtual void addPoints(const std::vector<std::vector<double> > & domain,
 //                         const std::vector<std::vector<double> > & response,
 //                         const std::vector<std::vector<std::vector<double> > > & grad)
 //     throw(SurfException);

   // add a point to the data set
   // 
   // point: an array of "pointCnt * dim" values representing the domain
   //        of the point to be added
   // fvals: an array of "pointCnt * respCnt" values representing the response
   //        values at this point
   // grads: an array of "pointCnt * dim * respCnt" values representing the
   //        gradient values at this point
   // dim: dimension of point (the dimension of the domain)
   // respCnt: the number of responses in this point
   // pointCnt: the number of points to be added
   //
   // 
 //  virtual void addPoints(const double * point,
 //                         const double * fvals, 
 //                         const double * grads,
 //                         int pointCnt,int dim,int respCnt);

   // add a set of data to this object
   //
   // data: SurfData object to have its data
   //       appended to this set of data
   //
 //  virtual void appendData(const SurfData & data);

   /// Return the dimensionality of the space from which the SurfPoints are drawn 
   virtual int getDimension() const;

   /// Return the number of response values in the data set
   virtual int getResponseCount() const;

   // get a vector to the points
   //
    virtual std::vector<SurfPoint*> & getPoints();

   /// Return an iterator for this class
   //void getIterator(AbstractDataIterator & iter) const;

   // get a iterator for this class
   //
   //void getIterator(AbstractDataIterator * iter) const;

   // so that you don't have to make the <<
   // operator a friend 
   //
   // os: stream to be written to
   //
   // *** what is a good format ***
   //virtual std::ostream & print(std::ostream & os) const;

   /// Read a set of SurfPoints from an input stream
   virtual std::istream& read(std::istream& is,bool binary=false);

   /// Write a set of SurfPoints to an output stream
   virtual std::ostream& write(std::ostream& os=cout, bool binary=false) const;

   /// Inform this object that a Surface wants to be notified
   /// whenever this object changes
   virtual void addListener(Surface *);
 
   /// remove the Surface from the list of surfaces that are notified
   // when the data changes
   virtual void removeListener(Surface *);

// overloaded operators
////////////////////////////////////
   
   // assignment operator
   //
   SurfData& operator=(const SurfData& sd);
    
};

// so a SurfData object can be printed
//
std::ostream & operator<<(std::ostream & os,const SurfData & data);

#endif
