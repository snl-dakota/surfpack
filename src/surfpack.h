#include "config.h"

#ifndef __SURFPACK_H__
#define __SURFPACK_H__

#include <cassert>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

//class AbstractSurfDataIterator;
class Surface;
class SurfData;

namespace surfpack {


// _____________________________________________________________________________
// Constants 
// _____________________________________________________________________________
  // Length of the field for double-precision number stream output
  const unsigned field_width = 26;

  // Precision of output for double precision numbers
  const unsigned output_precision = 17;

// _____________________________________________________________________________
// Nested Types 
// _____________________________________________________________________________

  /// For use in comparing the actual value of some function with the estimate
  /// given by a Surface approximation
  struct ErrorStruct {
    double observed;
    double estimated;
  };

  /// Thrown when an attempt to open a file for reading or writing fails
  class file_open_failure: public std::runtime_error
  {
  public:
    file_open_failure(const std::string& filename = "") 
      : std::runtime_error("File " + filename + " could not be opened.") {}
  };
    
  /// Thrown when end-of-file is reached unexpectedly, when an unrecognized or
  /// unacceptable file extension is encountered, or when a file contains
  /// unexpected or illegally formatted contents.
  class io_exception: public std::runtime_error
  {
  public:
    io_exception(const std::string& msg = "") : std::runtime_error(msg) {}
  };
  
// ____________________________________________________________________________
// I/O 
// ____________________________________________________________________________
  
  /// Write the value of contents to the file specified by filename.  Throw an
  /// exception if the file cannot be opened.
  void writeFile(std::string filename, std::string contents);

  /// Write the parameter header, followed by the matrix mat (the dimensions of
  /// which are specified by parameters rows and columns) to the parameter os.
  /// If c_style is true, the memory layout is assumed to follow the C
  /// convention (if mat points to an m by n matrix, the first m values are
  /// interpreted as the first row).  Otherwise, the layout is assumed to 
  /// follow the Fortran convention (the first n values are interpreted as the 
  /// first column).
  void writeMatrix(const std::string header, double* mat, unsigned rows, 
    unsigned columns, std::ostream& os, bool c_style = false);

  /// Write the parameter header, followed by the matrix mat (the dimensions of
  /// which are specified by parameters rows and columns) to the parameter os.
  /// If c_style is true, the memory layout is assumed to follow the C
  /// convention (if mat points to an m by n matrix, the first m values are
  /// interpreted as the first row).  Otherwise, the layout is assumed to 
  /// follow the Fortran convention (the first n values are interpreted as the 
  /// first column).
  void writeMatrix(const std::string header, unsigned* mat, unsigned rows, 
    unsigned columns, std::ostream& os, bool c_style = false);

  /// Write the contents of a matrix to a file specified by parameter filename.
  /// Open the file and call another version of writeMatrix.
  void writeMatrix(const std::string filename, double* mat, unsigned rows, 
    unsigned columns, bool c_style = false);

  /// Write the contents of a matrix to a file specified by parameter filename.
  /// Open the file and call another version of writeMatrix.
  void writeMatrix(const std::string filename, unsigned* mat, unsigned rows, 
    unsigned columns, bool c_style = false);

  /// Write the parameter header followed by the values in the vector
  void printVector(const std::string header, std::vector<double>& vec, 
    std::ostream& os = std::cout);

  /// Return true if the file specified by parameter file name has the extension
  /// specified by parameter extension
  bool hasExtension(const std::string& filename, const std::string extension);

  /// Throw an exception if end-of-file has been reached 
  void checkForEOF(std::istream& is);

  /// Open the file specified by filename and return the type of Surface that 
  /// is written there.  Throw an exception if the file cannot be opened, or if
  /// the file extension is anything other than .txt or .srf. 
  const std::string surfaceName(const std::string filename);

  /// Return the next item in the file as a string.  If the file is opened in
  /// binary mode, first read an integer that specifies the number of 
  /// characters in the string, then read the string. 
  const std::string readName(std::istream& is, bool binary);

// ____________________________________________________________________________
// Vector helper methods 
// ____________________________________________________________________________

  /// Return the arithmetic mean (average) of the values in vector vals
  double mean(std::vector<double>& vals);

  /// Return the sample variance of the values in vals
  double sample_var(std::vector<double>& vals);
  
  /// Return the sample standard deviation of the values in vals
  double sample_sd(std::vector<double>& vals);

  /// Return the sum of squared deviations from the mean
  double sum_squared_deviations(std::vector<double>& vals);
  
  /// Return the euclidean distance between pt1 and pt2.  Throw an exception if
  /// the dimensionality of the two vectors does not match.
  double euclideanDistance(const std::vector<double>& pt1, 
    const std::vector<double>& pt2);

  /// Store the vector difference between pt1 and pt2 in the paramter diff.
  /// Throw an exception if the dimensionality of the points does not match.
  void vectorDifference(std::vector<double>& diff, 
    const std::vector<double>& pt1, const std::vector<double>& pt2);

// ____________________________________________________________________________
// Testing 
// ____________________________________________________________________________
  /// Return the value of the test function specified by parameter name at the 
  /// point specified by parameter pt
  /// \todo Change if-else construct to lookup in STL map of <name, function>
  /// pairs.
  double testFunction(const std::string name, const std::vector<double>& pt);

  /// Non-trivial polynomial function
  double moderatepoly(const std::vector<double>& pt);

  /// Tony Giunta's test function: a little wave mixed with a big wave, plus 
  /// noise
  double quasisine(const std::vector<double>& pt);

  /// f(x) = sigma{i=1 to n}(x_i^2 - 10*cos(2*pi*x) + 10).  With side 
  /// constraints near zero (e.g. +/-10) along each dimension, the function 
  /// appears highly-multimodal.  With larger bounds, the bowl shape becomes
  /// more dominant and the waviness is reduced to noise.
  double rastrigin(const std::vector<double>& pt);

  /// A multi-dimensional extension of the classic Rosenbrock test function
  double rosenbrock(const std::vector<double>& pt);

  /// f(x) = 3 + sigma{i=1 to n}(2*x_i)
  double simplepoly(const std::vector<double>& pt);

  /// Sum of the sine function along each dimension
  double sinewave(const std::vector<double>& pt);

  /// Sum of squares along each dimension
  double sphere(const std::vector<double>& pt);

  /// f(x) = sigma{i=1 to i}(x_i)
  double sumofall(const std::vector<double>& pt);

  /// f(x) = sigma{i=1 to n}(x_i + sin x_i)
  double xplussinex(const std::vector<double>& pt);
} // namespace surfpack
#endif
