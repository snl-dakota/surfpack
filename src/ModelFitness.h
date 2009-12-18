/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#ifndef __MODEL_FITNESS_H__
#define __MODEL_FITNESS_H__

#ifdef HAVE_CONFIG_H
#include "surfpack_config.h"
#endif
#include "surfpack_system_headers.h"
#include "surfpack.h"

//class SurfPoint;
//class SurfData;
//class SurfScaler;
//namespace surfpack { class ErrorStruct; }
//

// BMA: For AIX included the surfpack.h above and commented enum
//enum MetricType;
template< typename T>
std::vector< T >& vecSubInPlace(std::vector< T >& sub, std::vector< T >& min)
{
  assert(sub.size() == min.size());
  typedef typename std::vector< T >::iterator VecIt;
  VecIt subIt, minIt;
  for (subIt = sub.begin(), minIt = min.begin();
	subIt != sub.end();
	++subIt, ++minIt) {
    *subIt -= *minIt;
  }
  return sub;
}

template< typename T>
std::vector< T > vecSub(std::vector< T >& sub, std::vector< T >& min)
{
  assert(sub.size() == min.size());
  std::vector< T > result(sub.size());
  typedef typename std::vector< T >::iterator VecIt;
  VecIt subIt, minIt, resIt;
  for (subIt = sub.begin(), minIt = min.begin(), resIt = result.begin();
	subIt != sub.end();
	++subIt, ++minIt, ++resIt) {
    *resIt = *subIt - *minIt;
  }
  return result;
}


class Residual
{
public:
  Residual(DifferenceType dt_in);
  double operator()(double observed, double predicted) const; 
protected:
  DifferenceType dt;
};

class VecSummary
{
public:
  VecSummary(MetricType mt_in);
  MetricType mt;
  double operator()(const VecDbl& resids) const;
};

class ModelFitness 
{
public:
  virtual double operator()(const SurfpackModel& sm, const SurfData& sd) const = 0;
  virtual double operator()(const VecDbl& obs, const VecDbl& pred) const;
  static ModelFitness* Create(const std::string& metric, unsigned n = 0);
  static VecDbl getResiduals(const Residual& resid, const VecDbl& obs, const VecDbl& pred);
};

class StandardFitness : public ModelFitness
{
public:
  StandardFitness();
  StandardFitness(const Residual& resid_in, const VecSummary& vecsumry_in);
  virtual double operator()(const SurfpackModel& sm, const SurfData& sd) const;
  virtual double operator()(const VecDbl& obs, const VecDbl& pred) const;
protected:
  Residual resid;
  VecSummary vecsumry;
};

class CrossValidationFitness : public ModelFitness
{
public:
  CrossValidationFitness(unsigned n_in);
  virtual double operator()(const SurfpackModel& sm, const SurfData& sd) const;
  unsigned n;
};

class PRESSFitness: public ModelFitness
{
public:
  PRESSFitness();
  virtual double operator()(const SurfpackModel& sm, const SurfData& sd) const;
};

class R2Fitness: public ModelFitness
{
public:
  R2Fitness();
  virtual double operator()(const SurfpackModel& sm, const SurfData& sd) const;
};

#endif
