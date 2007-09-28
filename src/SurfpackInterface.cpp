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
