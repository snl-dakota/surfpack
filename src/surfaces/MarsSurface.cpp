// These lines were automatically prepended
#ifndef SURFPACK_CONFIG_H
#define SURFPACK_CONFIG_H
#include "config.h"
#endif
// End of prepended lines

#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>
#include <set>
#include "surfpack.h"
#include "SurfPoint.h"
#include "SurfData.h"
#include "Surface.h"
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
using namespace surfpack;

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

MarsSurface::MarsSurface(SurfData* sd) : Surface(sd),
  xMatrix(0), fm(0), im(0)
{
#ifdef __TESTING_MODE__
  constructCount++;
#endif
}

MarsSurface::MarsSurface(const string filename) : Surface(0), xMatrix(0), 
  fm(0), im(0)
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
  if (sd) {
    return 4 * sd->xSize();
  } else {
    cerr << "Cannot compute minPointsRequired without data" << endl;
    return INT_MAX;
  }
}

double MarsSurface::evaluate(const std::vector<double>& x)
{
  //int nmcv = 0;
  //int ntcv = 0;
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

void MarsSurface::build(SurfData& data)
{
  delete [] xMatrix;
  delete [] fm;
  delete [] im;
  int nmcv = 0;
  int ntcv = 0;
  //int nk = 15;
  //int mi = 2;
  n = static_cast<int>(data.size());
  np = static_cast<int>(data.xSize());
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

  unsigned pts = data.size();
  for (unsigned i = 0; i < data.size(); i++) {
    //cout << current->getF(responseIndex) << endl;
    //const vector<double> domain = current.X();
    for (int j = 0; j < np; j++) {
      xMatrix[j*pts+i] = static_cast<real>(data[i].X()[j]); 
    }
    y[i] = static_cast<real>(data.getResponse(i));
    w[i] = 1.0f;
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
  originalData = true;
}

/// Create a surface of the same type as 'this.'  This objects data should
/// be replaced with the dataItr passed in, but all other attributes should
/// be the same (e.g., a second-order polynomial should return another 
/// second-order polynomial.  Surfaces returned by this method can be used
/// to compute the PRESS statistic.
MarsSurface* MarsSurface::makeSimilarWithNewData(SurfData* surfData)
{
  return new MarsSurface(surfData);
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
  np = static_cast<int>(sd->size());
  unsigned fmsize = 3+nk*(5*mi+nmcv+6)+2*np+ntcv;
  unsigned imsize = 21+nk*(3*mi+8);
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
    np = static_cast<int>(sd->size());
    unsigned fmsize = 3+nk*(5*mi+nmcv+6)+2*np+ntcv;
    unsigned imsize = 21+nk*(3*mi+8);
    
    unsigned old_precision = os.precision(surfpack::output_precision);
    os.setf(ios::scientific);
    os << setw(surfpack::field_width) << nmcv << "  nmcv" << endl
       << setw(surfpack::field_width) << ntcv << "  ntcv" << endl
       << setw(surfpack::field_width) << nk << "  nk" << endl
       << setw(surfpack::field_width) << mi << "  mi" << endl
       << setw(surfpack::field_width) << np << "  np" << endl;
    unsigned i;
    for (i = 0; i < fmsize; i++) {
      os << setw(surfpack::field_width) << fm[i] << "  fm[" << i << "]" << endl;
    }
    for (i = 0; i < imsize; i++) {
      os << setw(surfpack::field_width) << im[i] << "  im[" << i << "]" << endl;
    }
    os.unsetf(ios::scientific);
    os.precision(old_precision);
}

void MarsSurface::readBinary(std::istream& is)
{

  delete [] fm;
  delete [] im;
  int nmcv,ntcv,nk,mi,np;
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
  originalData = false;
}

void MarsSurface::readText(std::istream& is)
{
    delete [] fm;
    delete [] im;
    int nmcv,ntcv,nk,mi,np;
    string sline;
    istringstream streamline;
    getline(is,sline);
    streamline.str(sline); streamline.clear();
    streamline >> nmcv;
    getline(is,sline);
    streamline.str(sline); streamline.clear();
    streamline >> ntcv;
    getline(is,sline);
    streamline.str(sline); streamline.clear();
    streamline >> nk;
    getline(is,sline);
    streamline.str(sline); streamline.clear();
    streamline >> mi;
    getline(is,sline);
    streamline.str(sline); streamline.clear();
    streamline >> np;
    unsigned fmsize = 3+nk*(5*mi+nmcv+6)+2*np+ntcv;
    unsigned imsize = 21+nk*(3*mi+8);
    fm = new real[fmsize];
    im = new int[imsize];
    unsigned i;
    for (i = 0; i < fmsize; i++) {
      getline(is,sline);
      streamline.str(sline); streamline.clear();
      streamline >> fm[i];
    }
    for (i = 0; i < imsize; i++) {
      getline(is,sline);
      streamline.str(sline); streamline.clear();
      streamline >> im[i];
    }
    originalData = false;
}

//_____________________________________________________________________________
// Testing 
//_____________________________________________________________________________

#ifdef __TESTING_MODE__
  int MarsSurface::constructCount = 0;
  int MarsSurface::destructCount = 0;
#endif
