/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */
#include "surfpack_config.h"
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

/// Round values that are close to integers to integers
void surfpack::approximateByIntegers(vector<double>& vals, double epsilon)
{
  for(vector<double>::iterator iter = vals.begin(); iter != vals.end();
    ++iter) {
    double approx = static_cast<double>(static_cast<int>(*iter));
    if (abs(*iter-approx) < epsilon) {
      *iter = approx;
    }
  }
}
// ____________________________________________________________________________
// Vector helper methods 
// ____________________________________________________________________________

/// Return the sum of the vector of values
double surfpack::sum_vector(std::vector<double>& vals)
{
  double sum = 0;
  for (unsigned i = 0; i < vals.size(); i++) {
    sum += vals[i];
  }
  return sum;
}

/// Return the arithmetic mean (average) of the values in vector vals
double surfpack::mean(std::vector<double>& vals)
{
  return sum_vector(vals) / vals.size();
}

/// Return the sample variance of the values in vals
double surfpack::sample_var(std::vector<double>& vals)
{
  double sse = sum_squared_deviations(vals);
  return sse / (vals.size() - 1);
}

/// Return the sample standard deviation of the values in vals
double surfpack::sample_sd(std::vector<double>& vals)
{
  return sqrt(surfpack::sample_var(vals));
}

/// Return the sum of squared deviations from the mean
double surfpack::sum_squared_deviations(std::vector<double>& vals)
{
  double sse = 0;
  double avg = surfpack::mean(vals);
  for (unsigned i = 0; i < vals.size(); i++) {
    sse += (vals[i]-avg)*(vals[i]-avg);
  }
  return sse;
}
  
/// Return the sum of absolute deviations from the mean
double surfpack::sum_absolute_deviations(std::vector<double>& vals)
{
  double sae = 0;
  double avg = surfpack::mean(vals);
  for (unsigned i = 0; i < vals.size(); i++) {
    sae += (vals[i]-avg);
  }
  return sae;
}
/// Return absolute, squared, or relative differences of second and third
/// parameters through the first parameter
void surfpack::differences(std::vector<double>& results, 
  std::vector<double>& observed, std::vector<double>& predicted, 
  enum DifferenceType dp)
{
  results.resize(observed.size());
  for (unsigned i = 0; i < observed.size(); i++) {
    results[i] = abs(observed[i] - predicted[i]);
    switch (dp) {
      case SQUARED: results[i] *= results[i]; break;
      case SCALED: results[i] /= abs(observed[i]); break;
    }
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

// ____________________________________________________________________________
// Functions for common linear algebra tasks 
// ____________________________________________________________________________
void surfpack::linearSystemLeastSquares(SurfpackMatrix<double>& A, 
  std::vector<double>& x, std::vector<double>& b)
{
  // Rows in A must == size of b
  assert(A.getNRows() == b.size()); 
  // System must be square or over-constrained
  assert(A.getNRows() >= A.getNCols());
  int n_rows = static_cast<int>(A.getNRows());
  int n_cols = static_cast<int>(A.getNCols());
  // Client may supply a "blank" initialized vector for x
  int lwork = n_rows*n_cols * 2;
  vector<double> work(lwork);
  // values must be passed by reference to Fortran, so variables must be 
  // declared for info, nrhs, trans
  int info;
  int nrhs=1;
  char trans = 'N';
  DGELS_F77(trans,n_rows,n_cols,nrhs,&A[0][0],n_rows,&b[0],
    n_rows,&work[0],lwork,info);
  x = b;
  x.resize(n_cols);
}

void surfpack::leastSquaresWithEqualityConstraints(SurfpackMatrix<double>& A, 
  std::vector<double>& x, std::vector<double>& c,
  SurfpackMatrix<double>& B, std::vector<double>& d)
{
  int m = static_cast<int>(A.getNRows());
  int n = static_cast<int>(A.getNCols());
  int p = static_cast<int>(B.getNRows());
  int B_cols = static_cast<int>(B.getNCols());
  assert(B_cols == n);
  assert(p <= n);
  assert(n <= p + m);
  assert(x.size() == static_cast<unsigned>(n));
  int lwork = (m + n + p);
  lwork *= lwork;
  ///\todo Compute optimal blocksize before running dgglse
  vector<double> work(lwork);
  int info = 0;
  //SurfpackMatrix<double> Acopy = A;
  //SurfpackMatrix<double> Bcopy = B;
  //cout << endl;
  //cout << "A" << endl;
  //cout << A.asString() << endl;
  //cout << A.asArrayString() << endl;
  //cout << endl;
  //cout << "B" << endl;
  //cout << B.asString() << endl;
  //cout << B.asArrayString() << endl;
  //cout << "c" << endl;
  //copy(c.begin(),c.end(),ostream_iterator<double>(cout,"\n"));
  //cout << "d" << endl;
  //copy(d.begin(),d.end(),ostream_iterator<double>(cout,"\n"));
  //cout << "x before" << endl;
  //copy(x.begin(),x.end(),ostream_iterator<double>(cout,"\n"));
  DGGLSE_F77(m,n,p,&A[0][0],m,&B[0][0],p,&c[0],&d[0],&x[0],&work[0],lwork,info);
  //cout << "x after" << endl;
  //copy(x.begin(),x.end(),ostream_iterator<double>(cout,"\n"));
  //vector<double> result;
  //cout << "B*x after" << endl;
  //surfpack::matrixVectorMult(result,Bcopy,x);
  //copy(result.begin(),result.end(),ostream_iterator<double>(cout,"\n"));
  if (info != 0) throw string("Error in dgglse\n");
}

SurfpackMatrix< double >& surfpack::inverse(SurfpackMatrix< double>& matrix)
{
  int n_rows = static_cast<int>(matrix.getNRows());
  int n_cols = static_cast<int>(matrix.getNCols());
  /// \todo compute optimal blocksize
  int lwork = n_cols;  // should be optimal blocksize
  vector<int> ipvt(n_rows);
  vector<double> work(lwork);
  int lda = n_cols;
  int info = 0;
  DGETRF_F77(n_rows,n_cols,&matrix[0][0],lda,&ipvt[0],info);
  DGETRI_F77(n_rows,&matrix[0][0],lda,&ipvt[0],&work[0],lwork,info);
  return matrix;
}

SurfpackMatrix< double >& surfpack::LUFact(SurfpackMatrix< double>& matrix,
  vector<int>& ipvt)
{
  int n_rows = static_cast<int>(matrix.getNRows());
  int n_cols = static_cast<int>(matrix.getNCols());
  // dgetrf may reorder the rows, the mapping between old and new rows
  // is returned in ipvt.
  ipvt.resize(n_rows);
  int lda = n_cols;
  int info = 0;
  DGETRF_F77(n_rows,n_cols,&matrix[0][0],lda,&ipvt[0],info);
  return matrix;
}

SurfpackMatrix< double >& surfpack::inverseAfterLUFact(SurfpackMatrix< double>& matrix, vector<int>& ipvt)
{
  int n_rows = static_cast<int>(matrix.getNRows());
  int n_cols = static_cast<int>(matrix.getNCols());
  int lwork = n_cols;  // should be optimal blocksize
  vector<double> work(lwork);
  int lda = n_rows;
  int info = 0;
  DGETRI_F77(n_rows,&matrix[0][0],lda,&ipvt[0],&work[0],lwork,info);
  return matrix;
}

vector< double >& surfpack::matrixVectorMult(vector< double >& result,
  SurfpackMatrix< double >& matrix, vector< double >& the_vector)
{
  assert(matrix.getNCols() == the_vector.size());
  result.resize(matrix.getNRows());
  char transpose = 'N';
  int n_rows = static_cast<int>(matrix.getNRows());
  int n_cols = static_cast<int>(matrix.getNCols());
  int incx = 1;
  int incy = 1;
  double alpha = 1.0;
  double beta = 0.0;
  DGEMV_F77(transpose,n_rows,n_cols,alpha,&matrix[0][0],n_rows,&the_vector[0],
    incx, beta, &result[0], incy);
  return result;
}

double surfpack::dot_product(const std::vector< double >& vector_a, 
	     const std::vector< double >& vector_b)
{
  assert(vector_a.size() == vector_b.size()); 
  int size = static_cast<int>(vector_a.size());
  int inc = 1;
  // ddot will not violate the constness
  return DDOT_F77(size, const_cast<double*>( &vector_a[0] ), inc,
			const_cast<double*>( &vector_b[0] ), inc);
}

vector< double >& surfpack::vectorShift(vector< double >& the_vector,
  double shift_value)
{
  for (vector< double >::iterator iter = the_vector.begin();
       iter != the_vector.end(); ++iter) {
    *iter -= shift_value;
  }
  return the_vector;
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
  } else if (name == "noise") {
    return surfpack::noise(pt);
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
  for (unsigned i = 0; i < pt.size(); i++) {
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

/// Random (different queries for the same point will give different results)
double surfpack::noise(const std::vector<double>& pt)
{
  return static_cast<double>(rand());
}
