#include "config.h"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>
#include <set>
#include <cmath>
#include "surfpack.h"
#include "SurfPoint.h"
#include "SurfData.h"
#include "Surface.h"
#include "RBFNetSurface.h"

extern "C" void dgels_(char& trans, int& m, int& n, int& nrhs, double* A,
      int& lda, double* B, int& ldb, double* work, int& lwork, int& info);

extern "C" void dgemm_(char& transa, char& transb, int& m, int& n, int& k, 
  double& alpha, double* A, int& lda, double* B, int& ldb, double& beta, 
  double* C, int& ldc);

extern "C" void dgemv_(char& trans, int& m, int& n, double& alpha, 
  double* A, int& lda, double* x, int& incx, double& beta, double* y, 
  int& incy);

using namespace std;
//_____________________________________________________________________________
// Data members 
//_____________________________________________________________________________

const string RBFNetSurface::name = "RBFNet";

//_____________________________________________________________________________
// Creation, Destruction, Initialization
//_____________________________________________________________________________

RBFNetSurface::RBFNetSurface(SurfData* sd) : Surface(sd)
{
#ifdef __TESTING_MODE__
  constructCount++;
#endif
  init();
}

RBFNetSurface::RBFNetSurface(const string filename) : Surface(0)
{
#ifdef __TESTING_MODE__
  constructCount++;
#endif
  read(filename);
}
RBFNetSurface::~RBFNetSurface()
{
#ifdef __TESTING_MODE__
  destructCount++;
#endif
}

void RBFNetSurface::init()
{
  // .1 is a magic number.  More intelligent algorithm for choosing radii 
  // is needed 
  radius = 0.1;
}

//_____________________________________________________________________________
// Overloaded Operators 
//_____________________________________________________________________________

//_____________________________________________________________________________
// Queries
//_____________________________________________________________________________

const std::string RBFNetSurface::surfaceName() const
{
  return name;
}

unsigned RBFNetSurface::minPointsRequired() const
{
  if (sd) {
    return sd->xSize();
  } else {
    cerr << "Cannot determine min points required without data" << endl;
    return INT_MAX;
  }
}

double RBFNetSurface::evaluate(const std::vector<double>& x)
{
  //cout << "Evaluate----------------------------------" << endl;
  double result = 0.0;
  double temp;
  for (unsigned i = 0; i < centers.size(); i++) {
    temp = surfpack::euclideanDistance(x,centers[i].X());
    //cout << "dist: " << temp;
    // distance squared
    temp *= temp;
    // scale by size
    temp /= sizes[i];
    //cout << " center: " << i
    //     << " weight: " << weights[i]
    //     << " rbf: " << exp(-temp)
    //     << endl;
    result += weights[i]*exp(-temp);
  }
  //cout << "End Evaluate: " << result << "-------" << endl;
  return result;
}

//_____________________________________________________________________________
// Commands 
//_____________________________________________________________________________

void RBFNetSurface::build(SurfData& surfData)
{
  centers.clear();
  for (unsigned i = 0; i < surfData.size(); i++) {
    centers.push_back(SurfPoint(surfData[i].X()));
  }
  sizes.resize(surfData.size());
  for (unsigned j = 0; j < sizes.size(); j++) {
    sizes[j] = radius; 
      
  }
  int numpts = static_cast<int>(surfData.size());
  int numcenters = static_cast<int>(centers.size());
  //cout << "numpts: " << numpts
  //     << " numcenters: " << numcenters
  //     << endl;
  
  // allocate space for the matrices
  double* hMatrix = new double[numpts*numcenters];
  double* responseVector = new double[numpts];
  double* hTransHMatrix = new double[numcenters*numcenters];
  double* resultVector = new double[numcenters];

  // populate the H and y matrices 
  for (unsigned point = 0; point < surfData.size(); point++) {
    for(int centerIndex = 0; centerIndex < numcenters; centerIndex++) 
    {
      double temp = surfpack::euclideanDistance(surfData[point].X(),
        centers[centerIndex].X());
      temp *= temp;
      temp /= sizes[centerIndex];
      temp = exp(-temp);
      hMatrix[point+centerIndex*numpts] = temp; 
    }
    responseVector[point] = surfData.getResponse(point);
  }

  //writeMatrix("H Matrix",hMatrix,numpts,numcenters,cout);
  //writeMatrix("response Vector",responseVector,numpts,1,cout);
  // The equation that we want to solve is H'Hw = H'y (solving for w)
  // We just filled H and y
  // We first need to compute H'y
  
  // increment (in memory) for the vectors x & y in alpha*A*x + beta*y
  int inc = 1; 

  // specify transpose ('T') or no transpose ('N')
  // Since we want H'y, we need to do transpose
  char trans = 'T';

  // dgemv performs operation y <- alpha*A*x + beta*y (A can be transposed)
  // We have alpha=1 beta=0 (no scaling), so y <- A'x 
  // We want H'y, so for us it's y <- H'y; we need to put our y vector in both the x & y positions
  
  // no scaling
  double alpha = 1.0;

  // don't perform the optional addition
  double beta = 0.0; 

  dgemv_(trans,numpts,numcenters,alpha,hMatrix,numpts,responseVector,inc,
    beta,resultVector,inc);
  //writeMatrix("result vector after dgemv",resultVector,numcenters,1,cout);
  
  // now responseVector holds H'y

  // now we need to compute H'H
  // dgemm computes (ignoring the scalars) C <- A*B+C, where A is m x k, B is k x n, C is m x n
  // in our case A (hMatrix) is numpts x numcenters and hTransHMatrix is numcenters*numcenters
  int m = static_cast<int>(numcenters);
  int n = static_cast<int>(numcenters);
  int k = static_cast<int>(numpts);
  char transa = 'T';
  char transb = 'N';
  dgemm_(transa,transb,m,n,k,alpha,hMatrix,k,hMatrix,k,beta,hTransHMatrix,m);

  //writeMatrix("hTransHMatrix",hTransHMatrix,numcenters,numcenters,cout);
  // values must be passed by reference to Fortran, so variables must be declared for info, nrhs, trans
  int info;
  int nrhs=1;
  trans = 'N';
  //cout << "A Matrix: " << endl;
  ////writeMatrix(a,pts,coefficients.size(),cout);
  //cout << "B Vector: " << endl;
  ////writeMatrix(b,pts,1,cout);
  int lwork = numcenters * numcenters * 2;
  double* work = new double[lwork];
  dgels_(trans,numcenters,numcenters,nrhs,hTransHMatrix,numcenters,
    resultVector,numcenters,work,lwork,info);
  //writeMatrix("weights after dgels",resultVector,numcenters,1,cout);
  
  weights.resize(numcenters);
  for (int i = 0; i < numcenters; i++) {
    weights[i] = resultVector[i];
  }
  //cout << "A Matrix after: " << endl;
  ////writeMatrix(a,pts,numCoeff,cout);
  if (info < 0) {
          cerr << "dgels_ returned with an error" << endl;
  }
  delete [] work;
  delete [] resultVector;
  delete [] hMatrix;
  delete [] hTransHMatrix;
  delete [] responseVector;
}

void RBFNetSurface::config(const SurfpackParser::Arg& arg)
{
  string argname = arg.name;
  if (argname == "radius") {
    radius = arg.lval.real;
  } else {
    Surface::config(arg);
  }
}
/// Create a surface of the same type as 'this.'  This objects data should
/// be replaced with the dataItr passed in, but all other attributes should
/// be the same (e.g., a second-order polynomial should return another 
/// second-order polynomial.  Surfaces returned by this method can be used
/// to compute the PRESS statistic.
RBFNetSurface* RBFNetSurface::makeSimilarWithNewData(SurfData* sd)
{
  return new RBFNetSurface(sd);
}

//_____________________________________________________________________________
// Helper methods 
//_____________________________________________________________________________

//_____________________________________________________________________________
// I/O 
//_____________________________________________________________________________

void RBFNetSurface::writeBinary(std::ostream& os)
{
  unsigned numcenters = centers.size();
  unsigned xsize = centers[0].xSize();
  //unsigned fsize = centers[0].fSize();
  os.write(reinterpret_cast<char*>(&numcenters),sizeof(numcenters));
  os.write(reinterpret_cast<char*>(&xsize),sizeof(xsize));
  //os.write(reinterpret_cast<char*>(&fsize),sizeof(fsize));
  for (unsigned j = 0; j < numcenters; j++) {
    centers[j].writeBinary(os);
    os.write(reinterpret_cast<char*>(&sizes[j]),sizeof(sizes[j]));
    os.write(reinterpret_cast<char*>(&weights[j]),sizeof(weights[j]));
  }
}

void RBFNetSurface::writeText(std::ostream& os)
{
  ios::fmtflags old_flags = os.flags();
  unsigned old_precision = os.precision(surfpack::output_precision);
  os.setf(ios::scientific);
  os << centers.size() << " Number of centers" << endl
     << centers[0].xSize() << " Number of dimensions" << endl;
     //<< centers[0].fSize() << " Number of responses" << endl;
  for (unsigned j = 0; j < centers.size() ; j++) {
    centers[j].writeText(os);
    os << sizes[j] << " Radius of rbf #" << j << endl
       << weights[j] << " Weight of rbf #" << j << endl;
  }
  os.flags(old_flags);
  os.precision(old_precision);
}

void RBFNetSurface::readBinary(std::istream& is)
{
  unsigned numcenters;  
  unsigned xsize;
  //unsigned fsize;
  is.read(reinterpret_cast<char*>(&numcenters),sizeof(numcenters));
  is.read(reinterpret_cast<char*>(&xsize),sizeof(xsize));
  //is.read(reinterpret_cast<char*>(&fsize),sizeof(fsize));
  double size;
  double weight;
  centers.clear();
  for (unsigned j = 0; j < numcenters; j++) {
    centers.push_back(SurfPoint(xsize,0,is,true));
    is.read(reinterpret_cast<char*>(&size),sizeof(size));
    sizes.push_back(size);
    is.read(reinterpret_cast<char*>(&weight),sizeof(weight));
    weights.push_back(weight);
  }
}

void RBFNetSurface::readText(std::istream& is)
{
  unsigned numcenters;  
  unsigned xsize;
  // read numcenters
  string sline;
  getline(is,sline);
  istringstream streamline;
  streamline.str(sline);
  streamline >> numcenters;
  //cout << "numcenters: " << numcenters << endl;
  // read xsize (dimensionality of the rbf centers)
  getline(is,sline);
  streamline.str(sline);
  string buffer = streamline.str();
  streamline >> xsize;
  //cout << "xsize: " << xsize << endl;
  //cout << "buffer: " << buffer << endl;

  // read sizes and weights
  double size;
  double weight;
  for (unsigned j = 0; j < numcenters; j++) {
    // read the location of a center in as a surfpoint
    centers.push_back(SurfPoint(xsize,0,is,false));
    //centers[j].writeText(cout);
    //cout << endl;
    // read size of rbf center;
    getline(is,sline);
    streamline.str(sline);
    buffer = streamline.str();
    streamline >> size;
    //cout << "size: " << size << endl;
    //cout << "buffer: " << buffer << endl;
    
    sizes.push_back(size);
    // read weight of rbf center;
    getline(is,sline);
    streamline.str(sline);
    buffer = streamline.str();
    streamline >> weight;
    //cout << "weight: " << weight << endl;
    //cout << "buffer: " << buffer << endl;
    weights.push_back(weight);
  }
}

//_____________________________________________________________________________
// Testing 
//_____________________________________________________________________________

#ifdef __TESTING_MODE__
  int RBFNetSurface::constructCount = 0;
  int RBFNetSurface::destructCount = 0;
#endif
