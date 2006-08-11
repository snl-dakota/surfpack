/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#include "surfpack_config.h"
#include "surfpack.h"
#include "SurfData.h"
#include "MarsSurface.h"

using std::cerr;
using std::endl;
using std::ios;
using std::istream;
using std::istringstream;
using std::max;
using std::ostream;
using std::setw;
using std::string;
using std::vector;

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

const string MarsSurface::name = "mars";

//_____________________________________________________________________________
// Creation, Destruction, Initialization
//_____________________________________________________________________________

MarsSurface::MarsSurface(SurfData* sd) : Surface(sd),
  xMatrix(0), fm(0), im(0)
{
  init();
#ifdef __TESTING_MODE__
  constructCount++;
#endif
}

MarsSurface::MarsSurface(const string filename) : Surface(0), xMatrix(0), 
  fm(0), im(0)
{
  init();
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

void MarsSurface::init()
{
  max_bases = 25;
  max_interactions = 2;
  interpolation = 2;
}
//_____________________________________________________________________________
// Overloaded Operators 
//_____________________________________________________________________________

//_____________________________________________________________________________
// Queries
//_____________________________________________________________________________

const string MarsSurface::surfaceName() const
{
  return name;
}

unsigned MarsSurface::minPointsRequired() const
{
  if (xsize <= 0) {
    throw string(
      "Dimensionality of data needed to determine number of required samples."
    );
  } else {
    return xsize+1; 
  }
}

double MarsSurface::evaluate(const vector<double>& x)
{
  //int nmcv = 0;
  //int ntcv = 0;
  // inputs
  // m=1 for linear, 2 for cubic
  //int ngc = 1;
  //int ngs = 1;
  //int icx = 1;
  //int max_bases = 15;
  //int max_interactions = 2;
  //
  //// outputs from plot
  //// number of curves
  //int nc; 
  //double* crv = new double[ngc*2*max_bases];
  //// number of surfaces
  //int ns;
  //double* srf = new double[ngs*ngs*max_bases];
  //double* sp = new double[max(4*ngs*ngs,max(ngc,2*n))];
  //int* mm = new int[max(2*(max_interactions+1),nmcv)];
  //plot(m,this->x,fm,im,ngc,ngs,icx,nc,crv,ns,srf,sp,mm);
  
  //delete [] crv;
  //delete [] srf;
  //delete [] sp;
  //delete [] mm;
  int nval = 1;
  real* xVector = new real[static_cast<int>(x.size())];
  for (int i = 0; i < static_cast<int>(x.size()); i++) {
    xVector[i] = static_cast<real>(x[i]);
  }
  real* sp = new real[2];
  real* f = new real[1];
  FMODM_F77(interpolation,nval,xVector,fm,im,f,sp);
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
  //int max_bases = 15;
  //int max_interactions = 2;
  n = static_cast<int>(data.size());
  np = static_cast<int>(data.xSize());
  xMatrix = new real[n*np];
  real* y = new real[n];
  real* w = new real[n];
  int* lx = new int[np];
  fm = new real[3+max_bases*(5*max_interactions+nmcv+6)+2*np+ntcv];
  im = new int[21+max_bases*(3*max_interactions+8)];
  real* sp = new real[2*(n*(max(max_bases+1,2)+3)+
			 max(3*n+5*max_bases+np,max(2*np,4*n)))+
		      2*np+4*max_bases];
  double* dp = new double[2*(max(n*max_bases,(max_bases+1)*(max_bases+1))+
			     max((max_bases+2)*(nmcv+3),4*max_bases))];
  int* mm = new int[2*(n*np+2*max(max_interactions,nmcv))];

  //unsigned pts = data.size();
  //for (unsigned i = 0; i < pts; i++) {
  for (int i = 0; i < n; i++) {
    //cout << current->getF(responseIndex) << endl;
    //const vector<double> domain = current.X();
    for (int j = 0; j < np; j++) {
      xMatrix[j*n+i] = static_cast<real>(data[i].X()[j]); 
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
  MARS_F77(n,np,xMatrix,y,w,max_bases,max_interactions,lx,fm,im,sp,dp,mm);
  delete [] y;
  delete [] w;
  delete [] lx;
  delete [] sp;
  delete [] dp;
  delete [] mm;
}

void MarsSurface::config(const Arg& arg)
{
  string argname = arg.name;
  if (argname == "max_bases") {
    max_bases = arg.getRVal()->getInteger();
  } else if (argname == "max_interactions") {
    max_interactions = arg.getRVal()->getInteger();
  } else if (argname == "interpolation") {
    if (arg.getRVal()->getStringLiteral() == "linear") {
      interpolation = 1;
    } else if (arg.getRVal()->getStringLiteral() == "cubic") {
      interpolation = 2;
    } else {
      cerr << "Expected value for interpolation: 'linear' or 'cubic'" << endl;
    } 
  } else {
    Surface::config(arg);
  }
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

void MarsSurface::writeBinary(ostream& os)
{
  int nmcv = 0;
  int ntcv = 0;
  np = static_cast<int>(xsize);
  unsigned fmsize = 3+max_bases*(5*max_interactions+nmcv+6)+2*np+ntcv;
  unsigned imsize = 21+max_bases*(3*max_interactions+8);
  os.write(reinterpret_cast<char*>(&max_bases),sizeof(max_bases));
  os.write(reinterpret_cast<char*>(&max_interactions),sizeof(max_interactions));
  os.write(reinterpret_cast<char*>(&nmcv),sizeof(nmcv));
  os.write(reinterpret_cast<char*>(&ntcv),sizeof(ntcv));
  os.write(reinterpret_cast<char*>(&np),sizeof(np));
  os.write(reinterpret_cast<char*>(fm),fmsize*sizeof(fm[0]));
  os.write(reinterpret_cast<char*>(im),imsize*sizeof(im[0]));
}

void MarsSurface::writeText(ostream& os)
{
    int nmcv = 0;
    int ntcv = 0;
    np = static_cast<int>(xsize);
    unsigned fmsize = 3+max_bases*(5*max_interactions+nmcv+6)+2*np+ntcv;
    unsigned imsize = 21+max_bases*(3*max_interactions+8);
    
    unsigned old_precision = os.precision(surfpack::output_precision);
    os.setf(ios::scientific);
    os << setw(surfpack::field_width) << nmcv << "  nmcv" << endl
       << setw(surfpack::field_width) << ntcv << "  ntcv" << endl
       << setw(surfpack::field_width) << max_bases << "  max_bases" << endl
       << setw(surfpack::field_width) << max_interactions << "  max_interactions" << endl
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

void MarsSurface::readBinary(istream& is)
{

  delete [] fm;
  delete [] im;
  int nmcv,ntcv,npval;
  is.read(reinterpret_cast<char*>(&max_bases),sizeof(max_bases));
  is.read(reinterpret_cast<char*>(&max_interactions),sizeof(max_interactions));
  is.read(reinterpret_cast<char*>(&nmcv),sizeof(nmcv));
  is.read(reinterpret_cast<char*>(&ntcv),sizeof(ntcv));
  is.read(reinterpret_cast<char*>(&npval),sizeof(npval));
  unsigned fmsize = 3+max_bases*(5*max_interactions+nmcv+6)+2*npval+ntcv;
  unsigned imsize = 21+max_bases*(3*max_interactions+8);
  fm = new real[fmsize];
  im = new int[imsize];
  is.read(reinterpret_cast<char*>(fm),fmsize*sizeof(fm[0]));
  is.read(reinterpret_cast<char*>(im),imsize*sizeof(im[0]));
}

void MarsSurface::readText(istream& is)
{
    delete [] fm;
    delete [] im;
    int nmcv,ntcv,npval;
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
    streamline >> max_bases;
    getline(is,sline);
    streamline.str(sline); streamline.clear();
    streamline >> max_interactions;
    getline(is,sline);
    streamline.str(sline); streamline.clear();
    streamline >> npval;
    unsigned fmsize = 3+max_bases*(5*max_interactions+nmcv+6)+2*npval+ntcv;
    unsigned imsize = 21+max_bases*(3*max_interactions+8);
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
}

//_____________________________________________________________________________
// Testing 
//_____________________________________________________________________________

#ifdef __TESTING_MODE__
  int MarsSurface::constructCount = 0;
  int MarsSurface::destructCount = 0;
#endif
