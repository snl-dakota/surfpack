#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>
#include "SurfPoint.h"
#include "SurfData.h"
#include "Surface.h"
#include "SurfDataIterator.h"
#include "MarsSurface.h"

#ifndef RS6K
#define mars mars_ 
#define plot plot_ 
#define fmod fmod_
#endif

int MarsSurface::nk = 25;
int MarsSurface::mi = 2;

extern "C" void mars(int&, int&, real&, real&, real&, int&, int&, int&,
  real&, int&, real&, double&, int&);

extern "C" void fmod(int&, int&, real&, real&, int&, real&, real&);
using namespace std;

void printMatrix(real* mat, unsigned rows, unsigned columns, ostream& os)
{
  for (unsigned r = 0; r < rows; r++) {
    for (unsigned c = 0; c < columns; c++) {
      os << setw(15) << mat[r + c*rows];
    }
    os << endl;
  }
}

void printIntMatrix(int* mat, unsigned rows, unsigned columns, ostream& os)
{
  for (unsigned r = 0; r < rows; r++) {
    for (unsigned c = 0; c < columns; c++) {
      os << setw(15) << mat[r + c*rows];
    }
    os << endl;
  }
}
//_____________________________________________________________________________
// Data members 
//_____________________________________________________________________________

const string MarsSurface::name = "Mars";

//_____________________________________________________________________________
// Creation, Destruction, Initialization
//_____________________________________________________________________________

MarsSurface::MarsSurface(SurfData& sd, unsigned responseIndex) : Surface(0),
  xMatrix(0), fm(0), im(0)
{
#ifdef __TESTING_MODE__
  constructCount++;
#endif
  dataItr = new SurfDataIterator(sd, responseIndex);
  build();
}

MarsSurface::MarsSurface(AbstractSurfDataIterator* dataItr) : Surface(dataItr), xMatrix(0), fm(0),
  im(0)
{
#ifdef __TESTING_MODE__
  constructCount++;
#endif
  build();
}

MarsSurface::MarsSurface(const string filename) : Surface(0), xMatrix(0), fm(0), im(0)
{
#ifdef __TESTING_MODE__
  constructCount++;
#endif
  read(filename);
}
MarsSurface::~MarsSurface()
{
#ifdef __TESTING_MODE__
  destructCount++;
#endif
  delete [] xMatrix;
  delete [] fm;
  delete [] im;
}


//_____________________________________________________________________________
// Overloaded Operators 
//_____________________________________________________________________________

//_____________________________________________________________________________
// Queries
//_____________________________________________________________________________

const std::string MarsSurface::surfaceName() const
{
  return name;
}

unsigned MarsSurface::minPointsRequired() const
{
  return 4 * dataItr->xSize();
}

double MarsSurface::evaluate(const std::vector<double>& x)
{
  //int nmcv = 0;
  //int ntcv = 0;
  ensureValidity();
  // inputs
  // m=1 for linear, 2 for cubic
  int m = 2; 
  //int ngc = 1;
  //int ngs = 1;
  //int icx = 1;
  //int nk = 15;
  //int mi = 2;
  //
  //// outputs from plot
  //// number of curves
  //int nc; 
  //double* crv = new double[ngc*2*nk];
  //// number of surfaces
  //int ns;
  //double* srf = new double[ngs*ngs*nk];
  //double* sp = new double[max(4*ngs*ngs,max(ngc,2*n))];
  //int* mm = new int[max(2*(mi+1),nmcv)];
  //plot(m,this->x,fm,im,ngc,ngs,icx,nc,crv,ns,srf,sp,mm);
  
  //delete [] crv;
  //delete [] srf;
  //delete [] sp;
  //delete [] mm;
  int n = 1;
  real* xVector = new real[x.size()];
  for (unsigned i = 0; i < x.size(); i++) {
    xVector[i] = static_cast<real>(x[i]);
  }
  real* sp = new real[2];
  real* f = new real[1];
  fmod(m,n,xVector[0],fm[0],im[0],f[0],sp[0]);
  delete [] sp;
  delete [] xVector;
  real result = *f;
  delete [] f;
  return result;
}

//_____________________________________________________________________________
// Commands 
//_____________________________________________________________________________

void MarsSurface::build()
{
  if (!acceptableData()) {
    cerr << "Unacceptable data.  Could not build Mars Surface" << endl;
  } else {
    delete [] xMatrix;
    delete [] fm;
    delete [] im;
    int nmcv = 0;
    int ntcv = 0;
    //int nk = 15;
    //int mi = 2;
    n = static_cast<int>(dataItr->elementCount());
    np = static_cast<int>(dataItr->xSize());
    xMatrix = new real[n*np];
    real* y = new real[n];
    real* w = new real[n];
    int* lx = new int[np];
    fm = new real[ 3+nk*(5*mi+nmcv+6)+2*np+ntcv];
    im = new int[ 21+nk*(3*mi+8) ];
    real* sp = new real[2*(n*(max(nk+1,2)+3) + max(3*n+5*nk+np, max(2*np, 4*n))) 
      + 2*np + 4*nk];
    double* dp = new double[2*(max(n*nk,(nk+1)*(nk+1)) + max((nk+2)*(nmcv+3),4*nk))];
    int* mm = new int[2*(n*np+2*max(mi,nmcv))];

    dataItr->toFront();
    unsigned i = 0;
    unsigned pts = dataItr->elementCount();
    while (!dataItr->isEnd()) {
      SurfPoint& current = dataItr->currentElement();
      //cout << current->getF(responseIndex) << endl;
      //const vector<double> domain = current.X();
      for (int j = 0; j < np; j++) {
        xMatrix[j*pts+i] = static_cast<real>(current.X()[j]); 
      }
      y[i] = static_cast<real>(current.F(dataItr->responseIndex()));
      w[i] = 1.0f;
      dataItr->nextElement();
      i++;
    } 
    // Specify each variable to be 'unrestricted'
    for (int k = 0; k < np; k++) {
      lx[k] = 1;
    }
    //printMatrix(xMatrix,n,np,cout);
    //printMatrix(w,n,1,cout);
    //printMatrix(y,n,1,cout);
    //printIntMatrix(lx,np,1,cout);
    mars(n,np,xMatrix[0],y[0],w[0],nk,mi,lx[0],fm[0],im[0],sp[0],dp[0],mm[0]);
    delete [] y;
    delete [] w;
    delete [] lx;
    delete [] sp;
    delete [] dp;
    delete [] mm;
    valid = true;
    originalData = true;
  }
}

/// Create a surface of the same type as 'this.'  This objects data should
/// be replaced with the dataItr passed in, but all other attributes should
/// be the same (e.g., a second-order polynomial should return another 
/// second-order polynomial.  Surfaces returned by this method can be used
/// to compute the PRESS statistic.
MarsSurface* MarsSurface::makeSimilarWithNewData
  (AbstractSurfDataIterator* dataItr)
{
  return new MarsSurface(dataItr);
}

//_____________________________________________________________________________
// Helper methods 
//_____________________________________________________________________________

//_____________________________________________________________________________
// I/O 
//_____________________________________________________________________________

void MarsSurface::writeBinary(std::ostream& os)
{
  int nmcv = 0;
  int ntcv = 0;
  np = static_cast<int>(dataItr->xSize());
  //int nk = 15;
  //int mi = 2;
  unsigned fmsize = 3+nk*(5*mi+nmcv+6)+2*np+ntcv;
  unsigned imsize = 21+nk*(3*mi+8);
  unsigned nameSize = name.size();
  os.write(reinterpret_cast<char*>(&nameSize),sizeof(nameSize));
  os.write(name.c_str(),nameSize);
  os.write(reinterpret_cast<char*>(&nk),sizeof(nk));
  os.write(reinterpret_cast<char*>(&mi),sizeof(mi));
  os.write(reinterpret_cast<char*>(&nmcv),sizeof(nmcv));
  os.write(reinterpret_cast<char*>(&ntcv),sizeof(ntcv));
  os.write(reinterpret_cast<char*>(&np),sizeof(np));
  os.write(reinterpret_cast<char*>(fm),fmsize*sizeof(fm[0]));
  os.write(reinterpret_cast<char*>(im),imsize*sizeof(im[0]));
}

void MarsSurface::writeText(std::ostream& os)
{
    int nmcv = 0;
    int ntcv = 0;
    np = static_cast<int>(dataItr->xSize());
    //int nk = 15;
    //int mi = 2;
    unsigned fmsize = 3+nk*(5*mi+nmcv+6)+2*np+ntcv;
    unsigned imsize = 21+nk*(3*mi+8);

    os << nmcv << endl
       << ntcv << endl
       << nk << endl
       << mi << endl
       << np << endl;
    unsigned i;
    for (i = 0; i < fmsize; i++) {
      os << fm[i] << endl;
    }
    for (i = 0; i < imsize; i++) {
      os << im[i] << endl;
    }
}

void MarsSurface::readBinary(std::istream& is)
{

  delete [] fm;
  delete [] im;
  int nmcv,ntcv,nk,mi,np;
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
  is.read(reinterpret_cast<char*>(&nk),sizeof(nk));
  is.read(reinterpret_cast<char*>(&mi),sizeof(mi));
  is.read(reinterpret_cast<char*>(&nmcv),sizeof(nmcv));
  is.read(reinterpret_cast<char*>(&ntcv),sizeof(ntcv));
  is.read(reinterpret_cast<char*>(&np),sizeof(np));
  unsigned fmsize = 3+nk*(5*mi+nmcv+6)+2*np+ntcv;
  unsigned imsize = 21+nk*(3*mi+8);
  fm = new real[fmsize];
  im = new int[imsize];
  is.read(reinterpret_cast<char*>(fm),fmsize*sizeof(fm[0]));
  is.read(reinterpret_cast<char*>(im),imsize*sizeof(im[0]));
  valid = true;
  originalData = false;
}

void MarsSurface::readText(std::istream& is)
{
    delete [] fm;
    delete [] im;
    int nmcv,ntcv,nk,mi,np;
    is >> nmcv
       >> ntcv
       >> nk
       >> mi
       >> np;
    unsigned fmsize = 3+nk*(5*mi+nmcv+6)+2*np+ntcv;
    unsigned imsize = 21+nk*(3*mi+8);
    fm = new real[fmsize];
    im = new int[imsize];
    unsigned i;
    for (i = 0; i < fmsize; i++) {
      is >> fm[i];
    }
    for (i = 0; i < imsize; i++) {
      is >> im[i];
    }
    valid = true;
    originalData = false;
}

//_____________________________________________________________________________
// Testing 
//_____________________________________________________________________________

#ifdef __TESTING_MODE__
  int MarsSurface::constructCount = 0;
  int MarsSurface::destructCount = 0;
#endif
