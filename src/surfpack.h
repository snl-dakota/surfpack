#ifndef __SURFPACK_H__
#define __SURFPACK_H__

#include <string>
#include <vector>

struct ErrorStruct {
  double observed;
  double estimated;
};

class AbstractSurfDataIterator;
class Surface;
class SurfData;

const std::string surfaceName(const std::string filename);
Surface* createSurface(const std::string filename);
Surface* createSurface(const std::string type, 
  AbstractSurfDataIterator* dataItr);
Surface* createSurface(const std::string type,
  AbstractSurfDataIterator* dataItr, unsigned order);
Surface* createSurface(const std::string type, SurfData& sd, 
  unsigned responseIndex);
Surface* createSurface(const std::string type, SurfData& sd, 
  unsigned responseIndex, unsigned order);

double euclideanDistance(const std::vector<double>& pt1, const std::vector<double>& pt2);
void vectorDifference(std::vector<double>& diff, const std::vector<double>& pt1,
  const std::vector<double>& pt2);

void printVector(const std::string header, std::vector<double>& vec);

namespace surfpack {
  const unsigned field_width = 26;
  const unsigned output_precision = 17;
}

#endif
