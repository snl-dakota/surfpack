// These lines were automatically prepended
#ifndef SURFPACK_CONFIG_H
#define SURFPACK_CONFIG_H
#include "config.h"
#endif
// End of prepended lines

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "surfpack.h"

//______________________________________________________________________________
// Functions for use by Surface methods
//______________________________________________________________________________

using namespace std;

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

/// Write the parameter header followed by the values in the vector
/// \todo Use an output iterator instead.  Priority: very low.
void surfpack::printVector(const std::string header, vector<double>& vec)
{
  cout << header << " size: " << vec.size() << endl;
  for (unsigned i = 0; i < vec.size(); i++) {
    cout << i << " " << vec[i] << endl;
  }
}

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

/// Make sure eof has not been reached unexpectedly
void surfpack::checkForEOF(istream& is)
{
  if (is.eof()) {
    throw surfpack::io_exception("End of file reached unexpectedly.");
  }
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
    cout << header << endl;
  }
  for (unsigned r = 0; r < rows; r++) {
    for (unsigned c = 0; c < columns; c++) {
      if (c_style) {
        os << setw(15) << mat[c + r * columns];
      } else {
        os << setw(15) << mat[r + c * rows];
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
    cout << header << endl;
  }
  for (unsigned r = 0; r < rows; r++) {
    for (unsigned c = 0; c < columns; c++) {
      if (c_style) {
        os << setw(15) << mat[c + r * columns];
      } else {
        os << setw(15) << mat[r + c * rows];
      }
    }
    os << endl;
  }
}

void surfpack::writeMatrix(const string filename, double* mat, unsigned rows, 
  unsigned columns, bool c_style)
{
  ofstream outfile(filename.c_str(),ios::out);
  if (!outfile) {
    cerr << "Could not open file (" << filename << ") in writeMatrix." << endl;
    return;
  }
  writeMatrix("none",mat,rows,columns,outfile, c_style);
  outfile.close();
}

void surfpack::writeMatrix(const string filename, unsigned* mat, unsigned rows, 
  unsigned columns, bool c_style)
{
  ofstream outfile(filename.c_str(),ios::out);
  if (!outfile) {
    cerr << "Could not open file (" << filename << ") in writeMatrix." << endl;
    return;
  }
  writeMatrix("none",mat,rows,columns,outfile, c_style);
  outfile.close();
}

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

double surfpack::sphere(const vector<double>& pt) 
{
  double result = 0.0;
  for (unsigned i = 0; i < pt.size(); i++) {
    double x = pt[i];
    result += x*x;
  }
  return result;
}

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

double surfpack::sinewave(const std::vector<double>& pt)
{
  double result = 0.0;
  for (unsigned i = 0; i < pt.size(); i++) {
    double x = pt[i];
    result += sin(x);
  }
  return result;

}

double surfpack::xplussinex(const std::vector<double>& pt)
{
  double result = 0.0;
  for (unsigned i = 0; i < pt.size(); i++) {
    double x = pt[i];
    result += x + sin(x);
  }
  return result;

}

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

double surfpack::simplepoly(const std::vector<double>& pt)
{
  double result = 3.0;
  for (unsigned i = 0; i < pt.size(); i++) {
    double x = pt[i];
    result += 2.0*x;
  }
  return result;
}

double surfpack::rastrigin(const vector<double>& pt) 
{
  double result = 0.0;
  for (unsigned i = 0; i < pt.size(); i++) {
    double x = pt[i];
    result += x*x-10*cos(4.0*acos(0.0)*x)+10.0;
  }
  return result;
}

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

double surfpack::sumofall(const vector<double>& pt) 
{
  double result = 0.0;
  for (unsigned i = 0; i < pt.size() - 1; i++) {
    result += pt[i];
  }
  return result;
}

bool surfpack::hasExtension(const std::string& filename, const std::string extension)
{
  return (filename.find(extension) == filename.size() - extension.size());
}

