/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#ifdef HAVE_CONFIG_H
#include "surfpack_config.h"
#endif
#include "AxesBounds.h"
class SurfData;
class Surface;
class ParsedCommand;
class SurfpackParser;

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

