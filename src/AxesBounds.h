// These lines were automatically prepended
#ifndef SURFPACK_CONFIG_H
#define SURFPACK_CONFIG_H
#include "config.h"
#endif
// End of prepended lines

#ifndef POINT_DEFINITION_H
#define POINT_DEFINITION_H

class SurfData;
class SurfPoint;
#include <vector>
#include <string>


class AxesBounds
{

public:
  struct Axis {
      double min;
      double max;
      int pts;
      double interval;
  };

  AxesBounds( std::vector<Axis> axes_);
  AxesBounds( std::string filename);
  ~AxesBounds();
  SurfData* sampleGrid(const std::vector<std::string>& testFunctions);
  SurfData* sampleMonteCarlo(unsigned numPts, const std::vector<std::string>& testFunctions);
  //SurfData* populate();
  void initialize();
  void nextPoint();

protected:
  std::vector<int> point;
  std::vector<double> surfptx;
  std::vector<Axis> axes;
  unsigned ndims;
  unsigned npts;
};
  
#endif
