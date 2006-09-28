/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#ifndef __MARS_CPP_SURFACE_H__
#define __MARS_CPP_SURFACE_H__
#ifdef HAVE_CONFIG_H
#include "surfpack_config.h"
#endif
#include "Surface.h"

class Basis
{
public:
  Basis() {}
  virtual ~Basis() {}
  virtual double eval(const std::vector<double>& x) const = 0;
  virtual std::string asString() const = 0;
  virtual Basis* clone() const = 0;
};

class MarsBasis : public Basis
{
public:
  MarsBasis(double sign_in, double knot_in, unsigned var_index_in,
    const Basis* multiplier_in = 0);
  double eval(const std::vector<double>& x) const;
  std::string asString() const;
  Basis* clone() const;
private:
  double sign;
  double knot;
  unsigned var_index;
  const Basis* multiplier;
};
class UnityBasis : public Basis
{
public:
  UnityBasis();
  std::string asString() const;
  double eval(const std::vector<double>& x) const ;
  Basis* clone() const;
};

class MarsCppSurface : public Surface
{
//_____________________________________________________________________________
// Creation, Destruction, Initialization
//_____________________________________________________________________________

public:
  MarsCppSurface(SurfData* sd);
  MarsCppSurface(const std::string filename);
  ~MarsCppSurface();
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
  virtual MarsCppSurface* makeSimilarWithNewData(SurfData* surfData);

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
  std::vector<Basis*> bfs;
  std::vector<double> coeffs;
  static const std::string name;
  int n;
  int max_bases; 
  int max_interactions;
  int interpolation;
  

//_____________________________________________________________________________
// Testing 
//_____________________________________________________________________________

};

#endif
