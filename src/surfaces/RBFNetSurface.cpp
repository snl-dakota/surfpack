#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>
#include "surfpack.h"
#include "SurfPoint.h"
#include "SurfData.h"
#include "Surface.h"
#include "SurfDataIterator.h"
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

RBFNetSurface::RBFNetSurface(SurfData& sd, unsigned responseIndex) : Surface(0)
{
#ifdef __TESTING_MODE__
  constructCount++;
#endif
  dataItr = new SurfDataIterator(sd, responseIndex);
  init();
  build();
}

RBFNetSurface::RBFNetSurface(AbstractSurfDataIterator* dataItr) : Surface(dataItr)
{
#ifdef __TESTING_MODE__
  constructCount++;
#endif
  init();
  build();
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
  if (!dataItr) {
    cerr << "No data in RBFNetSurface::init" << endl;
    return;
  }
  for (unsigned i = 0; i < dataItr->elementCount(); i++) {
    centers.push_back(dataItr->currentElement());
    dataItr->nextElement();
  }
  sizes.resize(dataItr->elementCount());
  for (unsigned j = 0; j < sizes.size(); j++) {
    sizes[j] = 0.1;
  }
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
  return dataItr->xSize();
}

double RBFNetSurface::evaluate(const std::vector<double>& x)
{
  //cout << "Evaluate----------------------------------" << endl;
  ensureValidity();
  double result = 0.0;
  double temp;
  for (unsigned i = 0; i < centers.size(); i++) {
    temp = euclideanDistance(x,centers[i].X());
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

void RBFNetSurface::build()
{
  //cout << "RBFNetSurface::build----------------------------" << endl;
  if (!acceptableData()) {
    cerr << "Unacceptable data.  Could not build RBFNet Surface" << endl;
  } else {
    //unsigned xsize = dataItr->xSize();
    int numpts = static_cast<int>(dataItr->elementCount());
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
    dataItr->toFront();
    unsigned point = 0;
    while(!dataItr->isEnd()) {
      SurfPoint& sp = dataItr->currentElement();
      for(int centerIndex = 0; centerIndex < numcenters; centerIndex++) 
      {
        double temp = euclideanDistance(sp.X(),centers[centerIndex].X());
        temp *= temp;
        temp /= sizes[centerIndex];
        temp = exp(-temp);
        hMatrix[point+centerIndex*numpts] = temp; 
      }
      responseVector[point] = 
        dataItr->currentElement().F(dataItr->responseIndex());
      dataItr->nextElement();
      point++;
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
    valid = true;
    originalData = true;
    delete [] work;
    delete [] resultVector;
    delete [] hMatrix;
    delete [] hTransHMatrix;
    delete [] responseVector;
  }
}

/// Create a surface of the same type as 'this.'  This objects data should
/// be replaced with the dataItr passed in, but all other attributes should
/// be the same (e.g., a second-order polynomial should return another 
/// second-order polynomial.  Surfaces returned by this method can be used
/// to compute the PRESS statistic.
RBFNetSurface* RBFNetSurface::makeSimilarWithNewData
  (AbstractSurfDataIterator* dataItr)
{
  return new RBFNetSurface(dataItr);
}

//_____________________________________________________________________________
// Helper methods 
//_____________________________________________________________________________

//_____________________________________________________________________________
// I/O 
//_____________________________________________________________________________

void RBFNetSurface::writeBinary(std::ostream& os)
{
  unsigned nameSize = name.size();
  os.write(reinterpret_cast<char*>(&nameSize),sizeof(nameSize));
  os.write(name.c_str(),nameSize);
  unsigned numcenters = centers.size();
  unsigned xsize = centers[0].xSize();
  unsigned fsize = centers[0].fSize();
  os.write(reinterpret_cast<char*>(&numcenters),sizeof(numcenters));
  os.write(reinterpret_cast<char*>(&xsize),sizeof(xsize));
  os.write(reinterpret_cast<char*>(&fsize),sizeof(fsize));
  unsigned j;
  for (j = 0; j < numcenters; j++) {
    centers[j].writeBinary(os);
    os.write(reinterpret_cast<char*>(&sizes[j]),sizeof(sizes[j]));
    os.write(reinterpret_cast<char*>(&weights[j]),sizeof(weights[j]));
  }
}

void RBFNetSurface::writeText(std::ostream& os)
{

}

void RBFNetSurface::readBinary(std::istream& is)
{

  unsigned nameSize;
  is.read(reinterpret_cast<char*>(&nameSize),sizeof(nameSize));
  char* surfaceType = new char[nameSize+1];
  is.read(surfaceType,nameSize);
  surfaceType[nameSize] = '\0';
  string nameInFile(surfaceType);
  delete [] surfaceType;
  if (nameInFile != name) {
    cerr << "Surface name in file is not 'Kriging'." << endl;
    cerr << "Cannot build surface." << endl;
    return;
  }
  unsigned numcenters;  
  unsigned xsize;
  unsigned fsize;
  is.read(reinterpret_cast<char*>(&numcenters),sizeof(numcenters));
  is.read(reinterpret_cast<char*>(&xsize),sizeof(xsize));
  is.read(reinterpret_cast<char*>(&fsize),sizeof(fsize));
  //weights.resize(numcenters);
  //sizes.resize(numcenters);
  unsigned j;
  double size;
  double weight;
  for (j = 0; j < numcenters; j++) {
    centers.push_back(SurfPoint(xsize,fsize,is,true));
    is.read(reinterpret_cast<char*>(&size),sizeof(size));
    sizes.push_back(size);
    is.read(reinterpret_cast<char*>(&weight),sizeof(weight));
    weights.push_back(weight);
  }
  valid = true;
}

void RBFNetSurface::readText(std::istream& is)
{
    valid = true;
    originalData = false;
}

//_____________________________________________________________________________
// Testing 
//_____________________________________________________________________________

#ifdef __TESTING_MODE__
  int RBFNetSurface::constructCount = 0;
  int RBFNetSurface::destructCount = 0;
#endif
