/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#ifdef HAVE_CONFIG_H
#include "surfpack_config.h"
#endif

#include "surfpack.h"
#include "SurfData.h"
#include "Surface.h"
#include "SurfpackModel.h"
#include "ModelFitness.h"

using std::cout;
using std::endl;
using std::ifstream;
using std::ios;
using std::istringstream;
using std::max_element;
using std::min_element;
using std::ofstream;
using std::ostream;
using std::ostream_iterator;
using std::ostringstream;
using std::set;
using std::string;
using std::vector;

double squareSum(double& acc, double element) { return acc += element*element; }
double absSum(double& acc, double element) { return acc += fabs(element); }

StandardFitness::StandardFitness()
: resid(Residual(SQUARED)), vecsumry(MT_MEAN)
{

}

StandardFitness::StandardFitness(const Residual& resid_in, const VecSummary& vecsumry_in)
: resid(resid_in), vecsumry(vecsumry_in)
{

}

double StandardFitness::operator()(const SurfpackModel& sm, const SurfData& sd) const
{
  VecDbl observed = sm(sd);
  VecDbl predicted = sd.getResponses();
  VecDbl residuals = getResiduals(resid,observed,predicted);
  return vecsumry(residuals);
  //VecDbl residuals = vecSub(observed,predicted);
  //double sum = accumulate(residuals.begin(),residuals.end(),0.0,squareSum);
  //return sum;
}

ModelFitness* ModelFitness::Create(const std::string& metric)
{
  if (metric == "sum_squared") return new StandardFitness(Residual(SQUARED),VecSummary(MT_SUM));
  return new StandardFitness(Residual(SQUARED),VecSummary(MT_SUM));
}

Residual::Residual(DifferenceType dt_in) : dt(dt_in) 
{

}

double Residual::operator()(double observed, double predicted) const
{
    switch(dt) {
      case ABSOLUTE: return fabs(observed - predicted);
      case SCALED: return fabs(observed - predicted)/fabs(observed);
      case SQUARED: return (observed-predicted)*(observed-predicted);
    }
    assert(dt == ABSOLUTE || dt == SCALED || dt == SQUARED); 
    return 0.0;
}

VecSummary::VecSummary(MetricType mt_in) 
  : mt(mt_in) 
{

}

double VecSummary::operator()(const VecDbl& resids) const
{
  return surfpack::mean(resids);
}

VecDbl ModelFitness::getResiduals(const Residual& resid, const VecDbl& obs, const VecDbl& pred)
{
  assert(obs.size() == pred.size());
  VecDbl result(obs.size());
  for (unsigned i = 0; i < result.size(); i++) {
    result[i] = resid(obs[i],pred[i]);
  }
  return result;
}
