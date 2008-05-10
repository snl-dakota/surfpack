/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#ifndef __SURFPACK_MODEL_H__
#define __SURFPACK_MODEL_H__
#ifdef HAVE_CONFIG_H
#include "surfpack_config.h"
#endif

#include "surfpack_system_headers.h"
#include "ModelScaler.h"
#include "SurfpackMatrix.h"
#include "surfpack.h"

class SurfData;

///////////////////////////////////////////////////////////
///	Surfpack Model Parameters 
///////////////////////////////////////////////////////////

typedef std::pair< std::string, std::string > ModelParam;
typedef std::map< std::string, std::string> ParamMap;

///////////////////////////////////////////////////////////
///	Surfpack Model 
///////////////////////////////////////////////////////////

class SurfpackModel
{
public:
  SurfpackModel(unsigned ndims_in);
  SurfpackModel(const SurfpackModel& other);
  virtual VecDbl operator()(const SurfData& data) const;
  double operator()(const VecDbl& x) const;
  virtual VecDbl gradient(const VecDbl& x) const;
  virtual MtxDbl hessian(const VecDbl& x) const;
  virtual std::string asString() const = 0;
  /// return the value of some error metric
  double goodnessOfFit(const std::string MetricName, const SurfData& data); 
  ///double goodnessOfFit(const std::string MetricName); 
  /// Return R^2, which measures the proportion of variability in the data 
  /// accounted for by the model (the approximating surface).
  /// R^2 = 1-SSE/SST = SSR/SST. 
  double rSquared(const SurfData& data);
  /// For each point x in dataSet, construct a similar approximation Surface
  /// that includes all of the points in dataSet except x.  Then evaluate the
  /// Surface at x.  The difference between x and the estimate of x given the
  /// rest of the data is the residual for x.  The PRESS statistic is the
  /// square root of the mean of the squares of all the residuals.
  double press(const SurfData& data);
  double nFoldCrossValidation(const SurfData& data, unsigned n);
  /// Compute one of several goodness of fit metrics.  The observed parameter
  /// should be a list of observed (or true) function values; the vector of
  /// predicted values gives the corresponding estimates from this surface.
  /// The dt parameter specifies the kind of residuals to compute.  ABSOLUTE
  /// residuals are (observed - predicted), SQUARED residuals are the squares
  /// of the absolute residuals.  SCALED residuals are the ABSOLUTE residuals
  /// divided by the observed value.  Given the type of residuals, the client
  /// may request the min, max, sum, or mean of the set of residuals over all
  /// the given data points.  Two additional metrics are possible.  The
  /// relative maximum absolute error is the maximum absolute error divided
  /// by the standard deviation of the observed data.  The relative average
  /// absolute error is the mean absolute error divided by the standard
  /// deviation of observed data.
  double genericMetric(std::vector<double>& observed,
    std::vector<double>& predicted, enum MetricType mt, enum DifferenceType dt);
  virtual ~SurfpackModel();
  virtual void scaler(ModelScaler* ms);
  ModelScaler* scaler() const;
  unsigned size() const { return ndims;}
  const ParamMap& parameters() const { return args; }
  void parameters(const ParamMap& args) { this->args = args;}
protected:
  virtual double evaluate(const VecDbl& x) const = 0;
  unsigned ndims;
  ParamMap args;
  ModelScaler*  mScaler;
};

///////////////////////////////////////////////////////////
///	Surfpack Model Factory
///////////////////////////////////////////////////////////

class SurfpackModelFactory
{

public:
  SurfpackModelFactory();
  SurfpackModelFactory(const ParamMap& args);
  virtual SurfpackModel* Build(const SurfData& sd);
  virtual SurfpackModel* Create(const SurfData& sd) = 0;
  virtual SurfpackModel* Create(const std::string& model_string) = 0;
  /// the minimum number of points with which Surfpack will build a model
  virtual unsigned minPointsRequired();
  /// the recommended default number of points
  virtual unsigned recommendedNumPoints();
  virtual void config();
  const ParamMap& parameters() const;
  void add(const std::string& name, const std::string& value);
protected:
  ParamMap params;
  unsigned ndims;
  unsigned response_index;
};

#endif
