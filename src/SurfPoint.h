// ---------------------------------------------------------------------
// SurfPoint (old SurrogateBasePoint) is the base for points used by a Surface. 
// The point itself is stored as 'x'. We store the function value 'fvals' 
// at the point 
// ---------------------------------------------------------------------
//
// Project: SURFPACK++
//
// File:       SurfPoint.h
// Author:     Eric Cyr     
// Created:    ???
// Modified:   June 17, 2002
//             June 26-28, 2002
//             July 11, 2002
// Modified:   Mark Richards
// 	       June 24, 2003
//
// Description
// + SurfPoint class - a container class for a point, response value
//   and gradient information
// + Left shift (<<) operator for SurfPoint. Easy printing.
// ---------------------------------------------------------------------

#ifndef __SURF_POINT_H__
#define __SURF_POINT_H__

#include <iostream>
#include <iomanip>
#include <vector>

// INVARIANTS: for SurfPoint
// ---------------------------
// fvals.size() == grad.size() == responseCount
// x.size()==grad[i].size() || grad[i].size()==0 

/** A point in some real-valued domain space with one or more
 *  associated response variables. Gradient information may
 *  also be included. */
class SurfPoint {
  
protected:
	
  /// The point in the domain.  The size of the vector
  /// is the dimensionality of the space. 
  std::vector<double> x;          

  /// One or more function values at x (i.e., f1(x), f2(x) ... )
  std::vector<double> fvals;      

  /// Gradient information at x
//  std::vector<std::vector<double> > grad;

  /// Number of responses at each point
  int responseCount;              

public:

#ifdef __TESTING_MODE__
  friend class SurfPointUnitTest;

  static int constructCount;
  static int copyCount;
  static int destructCount;
#endif

  // Construtor. The point 'x' is set to a vector of zero-length, 'isf'
  // is set to false, and 'f' is set to zero.
 // SurfPoint();

  /// Read a point from an input stream.  The point has dimensionality nxvals, 
  /// with nfvals response variables and ngvals gradient values.
  SurfPoint(int nxvals, int nfvals, int ngvals, istream & is, bool binary=false);

  /// Copy Constructor. Copies 'x', 'isf', and 'f' from 'pt'.
  SurfPoint(const SurfPoint & pt);

  // Constructor with 'x' and 'f' specified. 
  // In this case, 'isf' is set to true.
  SurfPoint(const std::vector<double> & x_in, double f_in);

  //SurfPoint(const std::vector<double> & x_in, double f_in,const std::vector<double> &);

  SurfPoint(const std::vector<double> &,const std::vector<double> &);

  //SurfPoint(const std::vector<double> & x_in,const std::vector<double> & f_in,
  //          const std::vector<std::vector<double> > &);

  //SurfPoint(const double *,const double *,const double *,int,int);

  // Construtor with 'x' specified. 
  // In this case, 'isf' is set to false and 'f' to zero.
   SurfPoint(const std::vector<double> & x_in);

  // Virtual Destructor. 
  virtual ~SurfPoint();

  /// Return the dimensionality of the space 
  virtual int getDimension() const;

  /// Return the number of response variables
  virtual int getResponseCount() const;

  // Assignment operator. Copies 'x', 'isf' and 'f' from 'pt'.
  virtual SurfPoint& operator=(const SurfPoint& pt);

  /// write point to an output stream
  virtual void write(ostream & os, bool binary=false) const;

  /// comparison operator
  //virtual bool operator==(const SurfPoint& pt) const;

  // comparison operator
  //virtual bool operator<(const SurfPoint& pt) const;
 
  // comparison operator
  //virtual bool operator>(const SurfPoint& pt) const;

  /// Copy 'xvec' into 'x'. Does not affect 'f'.
  //virtual void setX(const std::vector<double>& xvec);

  /// Return a const reference to 'x'.
  virtual const std::vector<double> & getX() const;

  // Set 'isf' and 'f'. If 'isfval' is true, then 'isf' is set to true
  // and 'f' is set to 'fval'. Otherwise, 'isf' is set to false and
  // 'f' is set to zero. Does not affect 'x'.
  virtual void setF(int response, double fval);

  // set all the responses at once (deletes old responses)
  //virtual void setResponses(const std::vector<double> &);
  //virtual void setResponses(const std::vector<double> &,const std::vector<std::vector<double> > &);
  //virtual void setResponses(const double *,const double *,int respCnt);

  // add a single reponse
  //virtual void addResponse(double,const std::vector<double> &); 
  //virtual void addResponse(double,const double * =0); 
  virtual int addResponse(double val = 0); 

  // Get the function value of the point; that is, return 'f'.
  virtual double getF(int=0) const;

  // set the gradiant at this point
  //virtual void setGrad(const std::vector<double> & grad,int index=0);
 
  // get the gradiant at this point
  //virtual void getGrad(std::vector<double> & grad,int index=0) const;
  //virtual const std::vector<double> & getGrad(int index=0) const;

  // Output the contents of this object to the current stream with no
  // line feeds. Return a pointer to the stream. This is used in
  // conjunction with ostream to avoid making operator<< a friend of
  // this object.
  virtual std::ostream & print(std::ostream & stream) const;

};

std::ostream& operator<<(std::ostream& stream, const SurfPoint & pt); 

#endif
