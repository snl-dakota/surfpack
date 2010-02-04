/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#ifdef HAVE_CONFIG_H
#include "surfpack_config.h"
#endif
#include "surfpack_system_headers.h"
#include "SurfpackInterface.h"
#include "surfpack.h"
#include "AxesBounds.h"
#include "SurfpackParserArgs.h"
#ifndef DISABLE_STANDALONE_PARSER
#include "SurfpackInterpreter.h"
#endif
#include "SurfData.h"
#include "ModelFactory.h"

using std::cerr;
using std::cout;
using std::endl;
using std::ostream;
using std::runtime_error;
using std::string;
using std::vector;
using std::ostringstream;

#include "SurfpackModel.h"
#include "ModelFitness.h"
using SurfpackInterface::CreateAxes;
using SurfpackInterface::CreateSurface;
using SurfpackInterface::LoadData;
using SurfpackInterface::LoadModel;
using SurfpackInterface::Evaluate;

///////////////////////////////////////////////////////////////////////////////
/////		SurfpackInterface namespace functions			  /////
///////////////////////////////////////////////////////////////////////////////


SurfData* SurfpackInterface::LoadData(const std::string& filename)
{
  SurfData* sd = new SurfData(filename);
  assert(sd);
  return sd;
}

SurfData* SurfpackInterface::LoadData(const std::string& filename, unsigned n_predictors, unsigned n_responses, unsigned n_cols_to_skip)
{
  SurfData* sd = new SurfData(filename,n_predictors,n_responses,n_cols_to_skip);
  assert(sd);
  return sd;
}

SurfpackModel* SurfpackInterface::LoadModel(const std::string& filename)
{
  throw string("LoadModel is not currently supported");
  return 0;
}

void SurfpackInterface::Save(const SurfpackModel* model, const std::string& filename)
{
  ///\todo Add Write Method to Models
  //model->write(filename);
}

void SurfpackInterface::Save(const SurfData* data, const std::string& filename)
{
  data->write(filename);
}

SurfpackModel* SurfpackInterface::CreateSurface(const SurfData* sd, ParamMap& args)
{
  assert(sd);
  SurfpackModel* model = 0;
  // Surface* model = ModelFactory::Create(ParamMap)
  return model;
}

void SurfpackInterface::Evaluate(const SurfpackModel* model, SurfData* sd, 
  const std::string& response_name)
{
  assert(model);
  assert(sd);
  VecDbl responses = (*model)(*sd);
  sd->addResponse(responses, response_name);
}

void SurfpackInterface::Evaluate(SurfData* sd, const VecStr test_functions)
{
  assert(sd);
  for (VecStr::const_iterator itr = test_functions.begin();
    itr != test_functions.end(); ++itr) {
    VecDbl results(sd->size());
    for (unsigned i = 0; i < results.size(); i++) {
      results[i] = surfpack::testFunction(*itr,(*sd)(i));
    }
    sd->addResponse(results,*itr);
  } 
}

AxesBounds* SurfpackInterface::CreateAxes(const std::string axes)
{
  return new AxesBounds(axes);
}

SurfData* SurfpackInterface::CreateSample(const AxesBounds* axes, const VecUns grid_points)
{
  return axes->sampleGrid(grid_points);  
}

SurfData* SurfpackInterface::CreateSample(const AxesBounds* axes, unsigned n_samples)
{
  return axes->sampleMonteCarlo(n_samples);
}

double SurfpackInterface::Fitness(const SurfpackModel* model, SurfData* sd, 
const std::string& metric, unsigned response, unsigned n)
{
  assert(model);
  assert(sd);
  sd->setDefaultIndex(response);
  ModelFitness* mf = ModelFitness::Create(metric,n);
  double result = (*mf)(*model,*sd);
  delete mf;
  return result;
}

/// Doxygen comment
double SurfpackInterface::Fitness(const SurfpackModel*, const std::string& metric, 
unsigned response, unsigned n)
{
  throw string("Must pass data set to compute metric");
  return 0.0;
}


