#include "config.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "surfpack.h"

using namespace std;

// ____________________________________________________________________________
// I/O 
// ____________________________________________________________________________

/// Write the value of contents to the file specified by filename.  Throw an
/// exception if the file cannot be opened.
void surfpack::writeFile(std::string filename, std::string contents)
{
  ofstream outfile(filename.c_str(), ios::out);
  if (!outfile) {
    throw surfpack::file_open_failure(filename);
  }
  outfile << contents << endl;
  outfile.close();
}


/// Write the parameter header, followed by the matrix mat (the dimensions of
/// which are specified by parameters rows and columns) to the parameter os.
/// If c_style is true, the memory layout is assumed to follow the C
/// convention (if mat points to an m by n matrix, the first m values are
/// interpreted as the first row).  Otherwise, the layout is assumed to 
/// follow the Fortran convention (the first n values are interpreted as the 
/// first column).
void surfpack::writeMatrix(const string header, double* mat, unsigned rows, 
  unsigned columns, ostream& os, bool c_style)
{
  if (header != "none" && header != "") {
    os << header << endl;
  }
  for (unsigned r = 0; r < rows; r++) {
    for (unsigned c = 0; c < columns; c++) {
      if (c_style) {
        os << setw(surfpack::field_width) << mat[c + r * columns];
      } else {
        os << setw(surfpack::field_width) << mat[r + c * rows];
      }
    }
    os << endl;
  }
}

/// Write the parameter header, followed by the matrix mat (the dimensions of
/// which are specified by parameters rows and columns) to the parameter os.
/// If c_style is true, the memory layout is assumed to follow the C
/// convention (if mat points to an m by n matrix, the first m values are
/// interpreted as the first row).  Otherwise, the layout is assumed to 
/// follow the Fortran convention (the first n values are interpreted as the 
/// first column).
/// \todo Templatize the method.  Priority: low.
void surfpack::writeMatrix(const string header, unsigned* mat, unsigned rows, 
  unsigned columns, ostream& os, bool c_style)
{
  if (header != "none" && header != "") {
    os << header << endl;
  }
  for (unsigned r = 0; r < rows; r++) {
    for (unsigned c = 0; c < columns; c++) {
      if (c_style) {
        os << setw(surfpack::field_width) << mat[c + r * columns];
      } else {
        os << setw(surfpack::field_width) << mat[r + c * rows];
      }
    }
    os << endl;
  }
}

/// Write the contents of a matrix to a file specified by parameter filename.
/// Open the file and call another version of writeMatrix.
void surfpack::writeMatrix(const string filename, double* mat, unsigned rows, 
  unsigned columns, bool c_style)
{
  ofstream outfile(filename.c_str(),ios::out);
  if (!outfile) {
    throw file_open_failure(filename);
  }
  writeMatrix("none",mat,rows,columns,outfile, c_style);
  outfile.close();
}

/// Write the contents of a matrix to a file specified by parameter filename.
/// Open the file and call another version of writeMatrix.
void surfpack::writeMatrix(const string filename, unsigned* mat, unsigned rows, 
  unsigned columns, bool c_style)
{
  ofstream outfile(filename.c_str(),ios::out);
  if (!outfile) {
    throw file_open_failure(filename);
  }
  writeMatrix("none",mat,rows,columns,outfile, c_style);
  outfile.close();
}

/// Write the parameter header followed by the values in the vector
/// \todo Use an output iterator instead.  Priority: very low.
void surfpack::printVector(const std::string header, vector<double>& vec,
  ostream& os)
{
  os << header << " size: " << vec.size() << endl;
  for (unsigned i = 0; i < vec.size(); i++) {
    os << i << " " << vec[i] << endl;
  }
}

/// Return true if the file specified by parameter file name has the extension
/// specified by parameter extension
bool surfpack::hasExtension(const std::string& filename, const std::string extension)
{
  return (filename.find(extension) == filename.size() - extension.size());
}

/// Throw an exception if end-of-file has been reached 
void surfpack::checkForEOF(istream& is)
{
  if (is.eof()) {
    throw surfpack::io_exception("End of file reached unexpectedly.");
  }
}

/// Open the file specified by filename and return the type of Surface that 
/// is written there.  Throw an exception if the file cannot be opened, or if
/// the file extension is anything other than .txt or .srf. 
const string surfpack::surfaceName(const string filename)
{
  bool binary;
  if (surfpack::hasExtension(filename,".srf")) {
    binary = true;;
  } else if (surfpack::hasExtension(filename,".txt")) {
    binary = false;
  } else {
    throw surfpack::io_exception(
      "Unrecognized filename extension.  Use .srf or .txt"
    );
  }
  ifstream infile(filename.c_str(), (binary ? ios::in|ios::binary : ios::in));
  if (!infile) {
    throw surfpack::file_open_failure(filename);
  } 
  string nameInFile = readName(infile, binary);  
  infile.close();
  return nameInFile;
} 

/// Return the next item in the file as a string.  If the file is opened in
/// binary mode, first read an integer that specifies the number of 
/// characters in the string, then read the string. 
const string surfpack::readName(istream& is, bool binary)
{
  string nameInFile;
  if (!binary) {
    getline(is,nameInFile);
    return nameInFile;
  } else {
    unsigned nameSize;
    is.read(reinterpret_cast<char*>(&nameSize),sizeof(nameSize));
    char* surfaceType = new char[nameSize+1];
    is.read(surfaceType,nameSize);
    surfaceType[nameSize] = '\0';
    return string(surfaceType);
  }
}

// ____________________________________________________________________________
// Vector helper methods 
// ____________________________________________________________________________

/// Return the arithmetic mean (average) of the values in vector vals
double surfpack::mean(std::vector<double>& vals)
{
  double sum = 0;
  for (unsigned i = 0; i < vals.size(); i++) {
    sum += vals[i];
  }
  return static_cast<double>(sum) / vals.size();
}

/// Return the sample variance of the values in vals
double surfpack::sample_var(std::vector<double>& vals)
{
  double sse = 0;
  double avg = surfpack::mean(vals);
  for (unsigned i = 0; i < vals.size(); i++) {
    sse += (vals[i]-avg)*(vals[i]-avg);
  }
  return sse / (vals.size() - 1);
}

/// Return the sample standard deviation of the values in vals
double surfpack::sample_sd(std::vector<double>& vals)
{
  return sqrt(surfpack::sample_var(vals));
}

/// Return the euclidean distance between pt1 and pt2.  Throw an exception if
/// the dimensionality of the two vectors does not match.
double surfpack::euclideanDistance(const vector<double>& pt1, 
  const vector<double>& pt2)
{
  double distance = 0.0;
  if (pt1.size() != pt2.size()) {
    throw string(
      "Cannot compute euclidean distance.  Vectors have different sizes.");
  } else {
    for (unsigned i = 0; i < pt1.size(); i++) {
      distance += (pt1[i] - pt2[i]) * (pt1[i] - pt2[i]);
    }
    distance = sqrt(distance);
  }
  return distance;

}

/// Store the vector difference between pt1 and pt2 in the paramter diff.
/// Throw an exception if the dimensionality of the points does not match.
void surfpack::vectorDifference(vector<double>& diff, const vector<double>& pt1,
  const vector<double>& pt2)
{
  if (pt1.size() != pt2.size() || pt1.size() != diff.size()) {
    cerr << "Cannot compute vector difference: size mismatch." << endl;
    return;
  }
  for (unsigned i = 0; i < pt1.size(); i++) {
    diff[i] = pt1[i] - pt2[i];
  }
}

// ____________________________________________________________________________
// Testing 
// ____________________________________________________________________________

/// Return the value of the test function specified by parameter name at the 
/// point specified by parameter pt
/// \todo Change if-else construct to lookup in STL map of <name, function>
/// pairs.
double surfpack::testFunction(const string name, const vector<double>& pt)
{
  if (name == "rosenbrock") {
    return surfpack::rosenbrock(pt);
  } else if (name == "sphere") {
    return surfpack::sphere(pt);
  } else if (name == "sumofall") {
    return surfpack::sumofall(pt);
  } else if (name == "simplepoly") {
    return surfpack::simplepoly(pt);
  } else if (name == "moderatepoly") {
    return surfpack::moderatepoly(pt);
  } else if (name == "sinewave") {
    return surfpack::sinewave(pt);
  } else if (name == "quasisine") {
    return surfpack::quasisine(pt);
  } else if (name == "xplussinex") {
    return surfpack::xplussinex(pt);
  } else {
    return surfpack::rastrigin(pt);
  }
}

/// Non-trivial polynomial function
double surfpack::moderatepoly(const std::vector<double>& pt)
{
  double result = -3.0;
  for (unsigned i = 0; i < pt.size(); i++) {
    double x = pt[i];
    switch (i % 3) {
      case 0: result -= 2.0*(x-3.0); break;
      case 1: result += 1.0*(x+3.0)*(x+3.0); break;
      case 2: result += 2.0*(x-3.0)*(pt[(i+2)%3]); break;
    }
  }
  return result;
}

/// Tony Giunta's test function: a little wave mixed with a big wave, plus 
/// noise
double surfpack::quasisine(const std::vector<double>& pt)
{
  double result = 0.0;
  double c = 16.0/15.0;
  double e = 1.0;
  for (unsigned i = 0; i < pt.size(); i++) {
    double x = pt[i];
    result += sin(c*x-e) + sin(c*x-e)*sin(c*x-e) + .02*sin(40.0*(c*x-e));
  }
  return result;
}

/// f(x) = sigma{i=1 to n}(x_i^2 - 10*cos(2*pi*x) + 10).  With side 
/// constraints near zero (e.g. +/-10) along each dimension, the function 
/// appears highly-multimodal.  With larger bounds, the bowl shape becomes
/// more dominant and the waviness is reduced to noise.
double surfpack::rastrigin(const vector<double>& pt) 
{
  double result = 0.0;
  for (unsigned i = 0; i < pt.size(); i++) {
    double x = pt[i];
    result += x*x-10*cos(4.0*acos(0.0)*x)+10.0;
  }
  return result;
}

/// A multi-dimensional extension of the classic Rosenbrock test function
double surfpack::rosenbrock(const vector<double>& pt) 
{
  double result = 0.0;
  for (unsigned i = 0; i < pt.size() - 1; i++) {
    double x = pt[i];
    double xp = pt[i+1];
    result += 100.0*(xp-x*x)*(xp-x*x)+(x-1.0)*(x-1.0); 
  }
  return result;
}

/// f(x) = 3 + sigma{i=1 to n}(2*x_i)
double surfpack::simplepoly(const std::vector<double>& pt)
{
  double result = 3.0;
  for (unsigned i = 0; i < pt.size(); i++) {
    double x = pt[i];
    result += 2.0*x;
  }
  return result;
}

/// Sum of the sine function along each dimension
double surfpack::sinewave(const std::vector<double>& pt)
{
  double result = 0.0;
  for (unsigned i = 0; i < pt.size(); i++) {
    double x = pt[i];
    result += sin(x);
  }
  return result;
}

/// Sum of squares along each dimension
double surfpack::sphere(const vector<double>& pt) 
{
  double result = 0.0;
  for (unsigned i = 0; i < pt.size(); i++) {
    double x = pt[i];
    result += x*x;
  }
  return result;
}

/// f(x) = sigma{i=1 to i}(x_i)
double surfpack::sumofall(const vector<double>& pt) 
{
  double result = 0.0;
  for (unsigned i = 0; i < pt.size() - 1; i++) {
    result += pt[i];
  }
  return result;
}

/// f(x) = sigma{i=1 to n}(x_i + sin x_i)
double surfpack::xplussinex(const std::vector<double>& pt)
{
  double result = 0.0;
  for (unsigned i = 0; i < pt.size(); i++) {
    double x = pt[i];
    result += x + sin(x);
  }
  return result;
}
