// These lines were automatically prepended
#ifndef SURFPACK_CONFIG_H
#define SURFPACK_CONFIG_H
#include "config.h"
#endif
// End of prepended lines

//----------------------------------------------------------------------------
// Project: SURFPACK++
//
// File: 	ANNSurface.h
// Author: 	Mark Richards
//----------------------------------------------------------------------------

#ifndef __ANN_SURFACE_H__ 
#define __ANN_SURFACE_H__ 

class SurfData;
class Surface;
class ANNApprox;

typedef float real;

class ANNSurface : public Surface
{
//_____________________________________________________________________________
// Creation, Destruction, Initialization
//_____________________________________________________________________________

public:
  ANNSurface(SurfData* sd);
  ANNSurface(const std::string filename);
  ~ANNSurface();
  void init();

//_____________________________________________________________________________
// Overloaded Operators 
//_____________________________________________________________________________

//_____________________________________________________________________________
// Queries
//_____________________________________________________________________________

  virtual const std::string surfaceName() const;
  
  virtual unsigned minPointsRequired() const;
  
  virtual double evaluate(const std::vector<double>& x);
//_____________________________________________________________________________
// Commands 
//_____________________________________________________________________________

  virtual void build(SurfData& data);
  
  virtual void config(const SurfpackParser::ArgList& arglist);
  /// Create a surface of the same type as 'this.'  This objects data should
  /// be replaced with the dataItr passed in, but all other attributes should
  /// be the same (e.g., a second-order polynomial should return another 
  /// second-order polynomial.  Surfaces returned by this method can be used
  /// to compute the PRESS statistic.
  virtual ANNSurface* makeSimilarWithNewData(SurfData* sd);

//_____________________________________________________________________________
// Helper methods 
//_____________________________________________________________________________
  
//_____________________________________________________________________________
// I/O 
//_____________________________________________________________________________

  virtual void writeBinary(std::ostream& os);
  virtual void writeText(std::ostream& os);
  virtual void readBinary(std::istream& is);
  virtual void readText(std::istream& is);
//_____________________________________________________________________________
// Data members 
//_____________________________________________________________________________
protected:
  static const std::string name;
  ANNApprox* annObject;
  double norm_bound;
  double svdfactor;
  double percent;

//_____________________________________________________________________________
// Testing 
//_____________________________________________________________________________

#ifdef __TESTING_MODE__
public:
  static int constructCount;
  static int destructCount;
#endif

};

#endif
