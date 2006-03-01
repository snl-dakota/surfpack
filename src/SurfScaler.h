/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */
#ifndef SURF_SCALER_H
#define SURF_SCALER_H

#include "surfpack_config.h"
#include "surfpack_system_headers.h"

class SurfData;
#include "SurfPoint.h"
#include "SurfpackParserArgs.h"

class DimensionScaler
{
public:
  static const bool response_var = true;
  static const bool predictor_var = false;
  DimensionScaler();
  virtual ~DimensionScaler();
  virtual double scale(double value) = 0;
  virtual double descale(double value) = 0;
  virtual DimensionScaler* clone() const = 0;
  virtual std::string asString() = 0;
#ifdef __TESTING_MODE__
  friend class SurfScalerTest;
  friend class SurfDataTest;
  friend class SurfaceTest;
#endif
};

class NonScaler : public DimensionScaler
{
public:
  NonScaler();
  virtual ~NonScaler();
  virtual double scale(double value);
  virtual double descale(double value);
  virtual DimensionScaler* clone() const;
  virtual std::string asString();
};

class LogScaler : public DimensionScaler
{
public:
  LogScaler();
  virtual ~LogScaler();
  virtual double scale(double value);
  virtual double descale(double value);
  virtual DimensionScaler* clone() const;
  virtual std::string asString();
};

class NormalizingScaler : public DimensionScaler
{
public:
  NormalizingScaler(const SurfData& surf_data, unsigned dim, 
    bool response = DimensionScaler::predictor_var);
  NormalizingScaler(const NormalizingScaler& other);
  virtual ~NormalizingScaler();
  virtual double scale(double value);
  virtual double descale(double value);
  virtual DimensionScaler* clone() const;
  virtual std::string asString();
#ifdef __TESTING_MODE__
  friend class SurfScalerTest;
  friend class SurfDataTest;
  friend class SurfaceTest;
#endif
protected:
  /// To scale a value x, (compute x - offset) / divisor
  double offset;
  double divisor;
};

class SurfScaler
{
public:
// ____________________________________________________________________________
// Creation, Destruction, Initialization 
// ____________________________________________________________________________

  /// Initializes scaledPoint 
  SurfScaler();

  /// Makes deep copy
  SurfScaler(const SurfScaler& other);

  /// Delete members of scalers and responseScalers
  ~SurfScaler();
 
protected:
  /// Delete members of scalers and responseScalers
  void cleanup();

// ____________________________________________________________________________
// Queries 
// ____________________________________________________________________________

public:
  /// Return a scaled version of the parameter point 
  double scale(unsigned index, double value) const;

  ///// Return the scaled value of a response variable
  double scaleResponse(unsigned index, double value);

  ///// Invert the scaling on the response variable
  double descaleResponse(unsigned index, double value);

  /// Number of predictor variables that this object scales
  unsigned xSize();

  /// Number of responses that this object scales
  unsigned fSize();

  /// Express SurfScaler as a string, with one DimensionScaler per line
  std::string asString();

// ____________________________________________________________________________
// Commands 
// ____________________________________________________________________________

  /// Iterate over a data set, computing the scaling parameters for each 
  /// dimension.
  void normalizeAll(const SurfData& data_);

  /// Sets the scaler for a response or predictor variable.  
  void setDimensionScaler(unsigned index, const DimensionScaler& ds, 
    bool response = DimensionScaler::predictor_var);

  /// Use the same DimensionScaler for all variables
  void setAll(const DimensionScaler& ds);

  /// Let Surfpack guess how to best scale each variable
  void scaleAuto(const SurfData& sd);

  /// Use SurfpackParserArgs to determine scaling for each parameter
  void configList(const SurfData& sd, ArgList& args);
  
  /// Parse out info from arg to configure SurfScaler
  void config(const SurfData& sd, const Arg& arg);

  /// Size scalers and responseScalers to match data set; initialize
  /// each dim to non-scaling
  void sizeToMatch(const SurfData& sd);
  
// ____________________________________________________________________________
// Helper methods 
// ____________________________________________________________________________


protected:
  std::vector< DimensionScaler* > scalers;
  std::vector< DimensionScaler* > responseScalers;
   
#ifdef __TESTING_MODE__
  friend class SurfScalerTest;
  friend class SurfDataTest;
  friend class SurfaceTest;
#endif
};

#endif
