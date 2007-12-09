/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
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
#include "SurfpackInterpreter.h"
#include "SurfData.h"
#include "Surface.h"
#include "SurfaceFactory.h"

using std::cerr;
using std::cout;
using std::endl;
using std::ostream;
using std::runtime_error;
using std::string;
using std::vector;
using std::ostringstream;
using SurfpackInterface::CreateAxes;
using SurfpackInterface::CreateSurface;
using SurfpackInterface::Fitness;
using SurfpackInterface::Load;
using SurfpackInterface::Save;

// New Interface
#include "SurfpackModel.h"
#include "ModelFitness.h"
using NewInterface::CreateAxes;
using NewInterface::CreateSurface;
using NewInterface::LoadData;
using NewInterface::LoadModel;
using NewInterface::Evaluate;

///////////////////////////////////////////////////////////////////////////////
/////		SurfpackInterface namespace functions			  /////
///////////////////////////////////////////////////////////////////////////////

void SurfpackInterface::Load(SurfData*& data, const string filename)
{
  data = new SurfData(filename);
  assert(data);
}

void SurfpackInterface::Load(SurfData*& data, const string filename,
  unsigned n_vars, unsigned n_responses, unsigned skip_columns)
{
  data = new SurfData(filename,n_vars,n_responses,skip_columns);
  assert(data);
}

void SurfpackInterface::Load(Surface*& surface, const string filename)
{
  surface = SurfaceFactory::createSurface(filename);
}

void SurfpackInterface::Save(SurfData* data, const string filename)
{
  assert(data);
  data->write(filename);
}

void SurfpackInterface::Save(Surface* surface, const string filename)
{
  assert(surface);
  surface->write(filename);
}

void SurfpackInterface::CreateSurface(Surface*& surface, SurfData* data, const string type, int response_index)
{
  assert(data);
  data->setDefaultIndex(response_index);
  surface = SurfaceFactory::createSurface(type, data);
  surface->config(Arg::makeArg("xsize",data->xSize()));
}

void SurfpackInterface::Evaluate(Surface* surface, SurfData* data)
{
  surface->getValue(*data);
}

double SurfpackInterface::Fitness(Surface* surface, const string metric, 
 SurfData* data, int response_index)
{
  assert(surface);
  if (data) {
    data->setDefaultIndex(response_index);
  }
  return surface->goodnessOfFit(metric,data);
}

double SurfpackInterface::Fitness(Surface* surface, unsigned n, 
 SurfData* data, int response_index)
{
  assert(surface);
  if (data) {
    data->setDefaultIndex(response_index);
  }
  SurfData& valid_data = surface->checkData(data);
  return surface->nFoldCrossValidation(valid_data,n);
}

void SurfpackInterface::CreateAxes(AxesBounds*& ab, const string bounds)
{
  ab = new AxesBounds(bounds);
}

void SurfpackInterface::CreateSample(SurfData*& data, const AxesBounds& axes, 
  const vector<unsigned>& grid_points, const vector< string >& test_functions)
{
  data = axes.sampleGrid(grid_points,test_functions);
}

void SurfpackInterface::CreateSample(SurfData*& data, const AxesBounds& axes, 
    const std::vector<unsigned>& grid_points)
{
  vector<string> no_functions;
  data = axes.sampleGrid(grid_points,no_functions);
}

void SurfpackInterface::CreateSample(SurfData*& data, const AxesBounds& axes, 
	unsigned size, const vector< string >& test_functions)
{
  data = axes.sampleMonteCarlo(size, test_functions);
}

SurfData* SurfpackInterface::CreateSample(const std::string& bounds, const std::string& grid, const std::string& functions)
{
  AxesBounds ab(bounds);
  vector<unsigned> gridPts = surfpack::toVec<unsigned>(grid);
  vector<string> testFunctions = surfpack::toVec<string>(functions);
  return ab.sampleGrid(gridPts,testFunctions);
}

///////////////////////////////////////////////////////////////////////////////
/////		SurfpackInterface namespace functions			  /////
///////////////////////////////////////////////////////////////////////////////


SurfData* NewInterface::LoadData(const std::string& filename)
{
  SurfData* sd = new SurfData(filename);
  assert(sd);
  return sd;
}

SurfData* NewInterface::LoadData(const std::string& filename, unsigned n_predictors, unsigned n_responses, unsigned n_cols_to_skip)
{
  SurfData* sd = new SurfData(filename,n_predictors,n_responses,n_cols_to_skip);
  assert(sd);
  return sd;
}

SurfpackModel* NewInterface::LoadModel(const std::string& filename)
{

}

void NewInterface::Save(const SurfpackModel* model, const std::string& filename)
{
  ///\todo Add Write Method to Models
  //model->write(filename);
}

void NewInterface::Save(const SurfData* data, const std::string& filename)
{
  data->write(filename);
}

SurfpackModel* NewInterface::CreateSurface(const SurfData* sd, ParamMap& args)
{
  assert(sd);
  SurfpackModel* model = 0;
  // Surface* model = ModelFactory::Create(ParamMap)
  return model;
}

void NewInterface::Evaluate(const SurfpackModel* model, SurfData* sd, 
  const std::string& response_name)
{
  assert(model);
  assert(sd);
  VecDbl responses = (*model)(*sd);
  sd->addResponse(responses, response_name);
}

void NewInterface::Evaluate(SurfData* sd, const VecStr test_functions)
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

AxesBounds* NewInterface::CreateAxes(const std::string axes)
{
  return new AxesBounds(axes);
}

SurfData* NewInterface::CreateSample(const AxesBounds* axes, const VecUns grid_points)
{
  return axes->sampleGrid(grid_points);  
}

SurfData* NewInterface::CreateSample(const AxesBounds* axes, unsigned n_samples)
{
  return axes->sampleMonteCarlo(n_samples);
}

double NewInterface::Fitness(const SurfpackModel* model, SurfData* sd, 
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
double NewInterface::Fitness(const SurfpackModel*, const std::string& metric, 
unsigned response, unsigned n)
{
  throw string("Must pass data set to compute metric");
  return 0.0;
}


