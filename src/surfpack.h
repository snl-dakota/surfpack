#ifndef __SURFPACK_H__
#define __SURFPACK_H__

#include <string>
#include <vector>

struct ErrorStruct {
  double observed;
  double estimated;
};

//class AbstractSurfDataIterator;
class Surface;
class SurfData;

const std::string surfaceName(const std::string filename);
Surface* createSurface(const std::string& filename);
Surface* createSurface(const std::string& type, SurfData& surfData);
Surface* createSurface(const std::string& type, SurfData& surfData, unsigned order);

double euclideanDistance(const std::vector<double>& pt1, const std::vector<double>& pt2);
void vectorDifference(std::vector<double>& diff, const std::vector<double>& pt1,
  const std::vector<double>& pt2);

double mean(std::vector<double>& vals);
double sample_var(std::vector<double>& vals);
double sample_sd(std::vector<double>& vals);

void printVector(const std::string header, std::vector<double>& vec);

namespace surfpack {
  const unsigned field_width = 26;
  const unsigned output_precision = 17;
}

#endif
