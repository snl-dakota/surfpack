/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#ifndef __MARS_SURFACE_H__
#define __MARS_SURFACE_H__
#ifdef HAVE_CONFIG_H
#include "surfpack_config.h"
#endif
#include "Surface.h"

typedef float real;

#define MARS_F77 F77_FUNC(mars,MARS)
#ifdef __cplusplus
extern "C" /* prevent C++ name mangling */
#endif
void MARS_F77(int&, int&, real*, real*, real*, int&, int&, int*,
	      real*, int*, real*, double*, int*);

#define FMODM_F77 F77_FUNC(fmodm,FMODM)
#ifdef __cplusplus
extern "C" /* prevent C++ name mangling */
#endif
void FMODM_F77(int&, int&, real*, real*, int*, real*, real*);

class MarsSurface : public Surface
{
//_____________________________________________________________________________
// Creation, Destruction, Initialization
//_____________________________________________________________________________

public:
  MarsSurface(SurfData* sd);
  MarsSurface(const std::string filename);
  ~MarsSurface();
protected:
  void init();
public:

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
  
  virtual void config(const Arg& arg);
  /// Create a surface of the same type as 'this.'  This objects data should
  /// be replaced with the dataItr passed in, but all other attributes should
  /// be the same (e.g., a second-order polynomial should return another 
  /// second-order polynomial.  Surfaces returned by this method can be used
  /// to compute the PRESS statistic.
  virtual MarsSurface* makeSimilarWithNewData(SurfData* surfData);

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
  real* xMatrix;
  real* fm;
  int* im;
  int n;
  int np;
  int max_bases; 
  int max_interactions;
  int interpolation;

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
