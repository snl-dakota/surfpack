// These lines were automatically prepended
#ifndef SURFPACK_CONFIG_H
#define SURFPACK_CONFIG_H
#include "config.h"
#endif
// End of prepended lines

#include <iostream>
#include <string>
#include <vector>

class SurfData;

namespace SurfDataGenerator {

  typedef struct {
      double min;
      double max;
      int pts;
      double interval;
  } Axis;

  void readPointSpec(std::vector<Axis>& axes, std::istream& is);
  SurfData pointSpecToSurfData(std::string filename); 
  void pointSpecToSurfDataFile(std::string pointSpecFilename,
    std::string outputFilename); 
  void pointSpecToSurfDataFile(std::string pointSpecFilename,
    std::string outputFilename, std::string functionName); 
  void pointSpecToSurfDataFile(std::string pointSpecFilename,
    std::string outputFilename, std::vector<std::string> functionNames); 
  void printPointSpec(std::vector<Axis>& axes, std::ostream& os);

  class PointCounter {
    public:
      PointCounter(std::vector<Axis>& axes_);
      ~PointCounter();
      void initialize();
      void nextPoint();
      void print(std::ostream& os) const;
      bool atLastPoint() const;
      const double& operator[](unsigned i) const;
    private:
      std::vector<Axis>& axes;
      bool lastPoint;
      std::vector<double> point;
  };

} // namespace SurfDataGenerator
