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
#include "SurfpackModel.h"
#include "ModelFitness.h"
#include "ModelFactory.h"

using std::cout;
using std::endl;
using std::ifstream;
using std::ios;
using std::istringstream;
using std::max_element;
using std::min_element;
using std::random_shuffle;
using std::ofstream;
using std::ostream;
using std::ostream_iterator;
using std::ostringstream;
using std::set;
using std::string;
using std::vector;

double ModelFitness::operator()(const VecDbl& obs, const VecDbl& pred) const
{
  throw string("Not implemented for abstract ModelFitness class");
}

StandardFitness::StandardFitness()
: resid(Residual(SQUARED)), vecsumry(MT_MEAN)
{

}

StandardFitness::StandardFitness(const Residual& resid_in, const VecSummary& vecsumry_in)
: resid(resid_in), vecsumry(vecsumry_in)
{

}

double StandardFitness::operator()(const VecDbl& obs, const VecDbl& pred) const
{
  VecDbl residuals = getResiduals(resid,obs,pred);
  return vecsumry(residuals);
}

double StandardFitness::operator()(const SurfpackModel& sm, const SurfData& sd) const
{
  VecDbl predicted = sm(sd);
  VecDbl observed = sd.getResponses();
  VecDbl residuals = getResiduals(resid,observed,predicted);
  return vecsumry(residuals);
}

ModelFitness* ModelFitness::Create(const std::string& metric, unsigned n)
{
  if (metric == "sum_squared") {
    return new StandardFitness(Residual(SQUARED),VecSummary(MT_SUM));
  } else if (metric == "mean_squared") {
    return new StandardFitness(Residual(SQUARED),VecSummary(MT_MEAN));
  } else if (metric == "max_squared") {
    return new StandardFitness(Residual(SQUARED),VecSummary(MT_MAXIMUM));
  } else if (metric == "sum_scaled") {
    return new StandardFitness(Residual(SCALED),VecSummary(MT_SUM));
  } else if (metric == "mean_scaled") {
    return new StandardFitness(Residual(SCALED),VecSummary(MT_MEAN));
  } else if (metric == "max_scaled") {
    return new StandardFitness(Residual(SCALED),VecSummary(MT_MAXIMUM));
  } else if (metric == "sum_abs") {
    return new StandardFitness(Residual(ABSOLUTE),VecSummary(MT_SUM));
  } else if (metric == "mean_abs") {
    return new StandardFitness(Residual(ABSOLUTE),VecSummary(MT_MEAN));
  } else if (metric == "max_abs") {
    return new StandardFitness(Residual(ABSOLUTE),VecSummary(MT_MAXIMUM));
  } else if (metric == "cv") {
    return new CrossValidationFitness(n);
  } else if (metric == "rsquared") {
    return new R2Fitness();
  }
  string msg = "Metric " + metric + " not supported";
  throw msg; 
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
      default: assert(false);
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
  switch (mt) {
    case MT_SUM: return std::accumulate(resids.begin(),resids.end(),0.0);
    case MT_MEAN: return surfpack::mean(resids);
    case MT_MAXIMUM: 
      VecDbl::const_iterator itr = max_element(resids.begin(),resids.end());
      return *itr;
    //default: throw string("Unknown vec summary");
  }
  return 0.0;
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

CrossValidationFitness::CrossValidationFitness(unsigned n_in)
  : ModelFitness(), n(n_in)
{

}

double CrossValidationFitness::operator()(const SurfpackModel& sm, const SurfData& sd) const
{
  //cout << "CV Fitness: " << n << endl;
  SurfData my_data = sd; // Get non const copy
  ParamMap args = sm.parameters();
  VecUns indices(my_data.size()); 
  for (unsigned i = 0; i < indices.size(); i++) indices[i] = i;
  random_shuffle(indices.begin(),indices.end());
  VecDbl estimates(my_data.size());
  for (unsigned partition = 0; partition < n; partition++) {
    //cout << "part: " << partition << endl;
    SetUns excludedPoints;
    unsigned low = surfpack::block_low(partition,n,my_data.size());
    unsigned high = surfpack::block_high(partition,n,my_data.size());
    //cout << "low/high: " << low << " " << high << endl;
    for (unsigned k = low; k <= high; k++) excludedPoints.insert(indices[k]);
    my_data.setExcludedPoints(excludedPoints);
    //cout << " excludes: " << excludedPoints.size() << endl;
    SurfpackModelFactory* factory = ModelFactory::createModelFactory(args);
    assert(my_data.size() >= factory->minPointsRequired());
    SurfpackModel* model = factory->Build(my_data);
    my_data.setExcludedPoints(SetUns());
    for (unsigned k = low; k <= high; k++) {
      estimates[indices[k]] = (*model)(my_data(k));
      //cout << "k: " << estimates[k] << endl;
    }
    delete model;
    delete factory;
  }
  ModelFitness* mf = ModelFitness::Create("mean_squared");
  VecDbl responses = sd.getResponses();
  assert(responses.size() == estimates.size());
  double fitness = (*mf)(estimates,responses);
  //cout << "CV vals: " << surfpack::fromVec<double>(estimates) << endl;
  delete mf;
  return fitness;
}

R2Fitness::R2Fitness()
{

}

double R2Fitness::operator()(const SurfpackModel& sm, const SurfData& sd) const
{

  VecDbl predicted = sm(sd);
  VecDbl observed = sd.getResponses();
  double obs_mean = surfpack::mean(observed);
  VecDbl vec_mean = VecDbl(observed.size(),obs_mean);
  StandardFitness sum_squares = StandardFitness(Residual(SQUARED),VecSummary(MT_SUM));
  return sum_squares(predicted,vec_mean)/sum_squares(observed,vec_mean);
}

