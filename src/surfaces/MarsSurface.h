//----------------------------------------------------------------------------
// Project: SURFPACK++
//
// File: 	MarsSurface.h
// Author: 	Mark Richards
//----------------------------------------------------------------------------

#ifndef __MARS_SURFACE_H__
#define __MARS_SURFACE_H__

class SurfData;
class Surface;
class AbstractSurfDataIterator;

class MarsSurface : public Surface
{
//_____________________________________________________________________________
// Creation, Destruction, Initialization
//_____________________________________________________________________________

public:
  MarsSurface(SurfData& sd, unsigned responseIndex = 0);
  MarsSurface(AbstractSurfDataIterator* dataItr);
  MarsSurface(const std::string filename);
  ~MarsSurface();

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

  virtual void build();
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
  float* xMatrix;
  float* fm;
  int* im;
  int n;
  int np;

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
