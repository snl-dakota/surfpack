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
#include "AxesBounds.h"
class SurfData;
class Surface;
class ParsedCommand;
class SurfpackParser;
class SurfpackModel;

namespace SurfpackInterface
{
  void Load(SurfData*& data, const std::string filename);
  void Load(SurfData*& data, const std::string filename,
    unsigned n_vars, unsigned n_responses, unsigned skip_columns);
  void Load(Surface*& surface, const std::string filename);
  void Save(SurfData* data, const std::string filename);
  void Save(Surface* surface, const std::string filename);
  void CreateSurface(Surface*& surface, SurfData* data, const std::string type, 
    int response_index = 0);
  void Evaluate(Surface* surface, SurfData* data);
  double Fitness(Surface* surface, const std::string metric, 
    SurfData* data = 0, int response_index = 0);
  double Fitness(Surface* surface, unsigned n, 
    SurfData* data = 0, int response_index = 0);
  void CreateAxes(AxesBounds*&, const std::string data);
  void CreateSample(SurfData*& data, const AxesBounds& axes, 
    const std::vector<unsigned>& grid_points, 
    const std::vector<std::string>& test_functions);
  SurfData* CreateSample(const std::string& bounds, const std::string& grid, const std::string& functions);
  void CreateSample(SurfData*& data, const AxesBounds& axes, 
    const std::vector<unsigned>& grid_points); 
  void CreateSample(SurfData*& data, const AxesBounds& axes, 
    unsigned size, const std::vector<std::string>& test_functions);
};

namespace NewInterface
{
  SurfData* LoadData(const std::string& filename);
  SurfData* LoadData(const std::string& filename, unsigned n_predictors,
    unsigned n_responses, unsigned n_cols_to_skip);
  SurfpackModel* LoadModel(const std::string& filename);
  void Save(const SurfpackModel* model, const std::string& filename);
  void Save(const SurfData* data, const std::string& filename);
  SurfpackModel* CreateSurface(const SurfData* sd, ParamMap& args);
  void Evaluate(const SurfpackModel* model, SurfData* sd, 
    const std::string& response_name = "");
  void Evaluate(SurfData* sd, const VecStr test_functions);
  AxesBounds* CreateAxes(const std::string axes);
  SurfData* CreateSample(const AxesBounds* axes, const VecUns grid_points);
  SurfData* CreateSample(const AxesBounds* axes, unsigned n_samples);
  double Fitness(const SurfpackModel*, SurfData* sd, 
    const std::string& metric, unsigned response = 0, unsigned n = 0);
  double Fitness(const SurfpackModel*, const std::string& metric, 
    unsigned response = 0, unsigned n = 0);
};


// SurfData* LoadData(string filename)
// SurfData* LoadData(string filename, vars, resp, skips)

// SurfpackModel* LoadModel(string filename)

// void Save(SurfpackModel*, string filename)
// void Save(SurfData*, string filename)

// SurfpackModel* Create(SurfData& sd, ArgList& al) 

// void Evaluate(SurfpackModel*, SurfData*)

// AxesBounds* CreateAxes(string axes)

// SurfData* CreateSample(AxesBounds* axes, grid_points) 
// SurfData* CreateSample(Axesbounds* axes, num_samples) 

// double Fitness(SurfpackModel*, data, metric, response = 0)
// double Fitness(SurfpackModel*, metric, response = 0)

 
