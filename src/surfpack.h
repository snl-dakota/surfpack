// These lines were automatically prepended
#ifndef SURFPACK_CONFIG_H
#define SURFPACK_CONFIG_H
#include "config.h"
#endif
// End of prepended lines

#ifndef __SURFPACK_H__
#define __SURFPACK_H__

#include <string>
#include <vector>
#include <stdexcept>

struct ErrorStruct {
  double observed;
  double estimated;
};

//class AbstractSurfDataIterator;
class Surface;
class SurfData;


double euclideanDistance(const std::vector<double>& pt1, const std::vector<double>& pt2);
void vectorDifference(std::vector<double>& diff, const std::vector<double>& pt1,
  const std::vector<double>& pt2);

double mean(std::vector<double>& vals);
double sample_var(std::vector<double>& vals);
double sample_sd(std::vector<double>& vals);

void printVector(const std::string header, std::vector<double>& vec);

namespace surfpack {
  const std::string surfaceName(const std::string filename);
  const std::string readName(std::istream& is, bool binary);
  const unsigned field_width = 26;
  const unsigned output_precision = 17;

  class file_open_failure: public std::runtime_error
  {
  public:
    file_open_failure(const std::string& filename = "") 
      : std::runtime_error("File " + filename + " could not be opened.") {}
  };
    
  class io_exception: public std::runtime_error
  {
  public:
    io_exception(const std::string& msg = "") : std::runtime_error(msg) {}
  };
  
  /// Make sure eof has not been reached unexpectedly
  void checkForEOF(std::istream& is);
  
  void writeMatrix(const std::string header, double* mat, unsigned rows, 
    unsigned columns, std::ostream& os, bool c_style = false);
  void writeMatrix(const std::string header, unsigned* mat, unsigned rows, 
    unsigned columns, std::ostream& os, bool c_style = false);
  void writeMatrix(const std::string filename, double* mat, unsigned rows, 
    unsigned columns, bool c_style = false);
  void writeMatrix(const std::string filename, unsigned* mat, unsigned rows, 
    unsigned columns, bool c_style = false);

  double rosenbrock(const std::vector<double>& pt);
  double rastrigin(const std::vector<double>& pt);
  double sphere(const std::vector<double>& pt);
  double sumofall(const std::vector<double>& pt);
  bool hasExtension(const std::string& filename, const std::string extension);
} // namespace surfpack

#endif
