/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */
#ifndef __POLYNOMIAL_SURFACE_H__
#define __POLYNOMIAL_SURFACE_H__
#include "surfpack_config.h"
#include "SurfpackMatrix.h"
#include "Surface.h"

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
  virtual void config(const Arg& arg);
  
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

  PolynomialSurface& operator=(const PolynomialSurface& other);
  virtual void build(SurfData& data);


  /// Set the degree of the polynomial fit (e.g., linear=1, quadratic=2, etc.)
  void setOrder(unsigned order_in);


// ____________________________________________________________________________
// Helper methods 
// ____________________________________________________________________________

  //static unsigned fact(unsigned x);
  //static unsigned nChooseR(unsigned n, unsigned r); 
  void resetTermCounter() const;
  double computeTerm(const std::vector<double>& x) const;
  void accumulateLikeFactors(std::vector<unsigned>& factorCounts);
  void nextTerm() const;

  double computeDerivTerm(const std::vector<double>& x,
  const std::vector<unsigned>& factorCounts, 
  const std::vector<unsigned>& differentiationCounts) const;

  void gradient(const std::vector<double> & x, 
    std::vector<double>& gradient_vector);

  void hessian(const std::vector<double> & x, 
    SurfpackMatrix<double>& hessian);

  void setEqualityConstraints(unsigned asv,const SurfPoint& sp, double valuePtr,
    std::vector<double>* gradientPtr, SurfpackMatrix<double>* hessianPtr);
// ____________________________________________________________________________
// Data members 
// ____________________________________________________________________________

protected:
  static const std::string name;
  unsigned order;
  std::vector<double> coefficients;
  mutable std::vector<unsigned> digits;
  //std::vector<unsigned> factorCounts;
  SurfpackMatrix<double> eqConLHS;
  std::vector<double> eqConRHS;
private:
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
