#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>
#include "SurfData.h"
#include "SurfDataIterator.h"
#include "PolynomialSurface.h"

extern "C" void dgels_(char&, int&, int&, int&, double*,
                       int&, double*, int&, double*, int&, int&);

extern "C" double ddot_(int&, double*, int&, double*, int&);

extern "C" void dgemv_(char&, int&, int&, double&, double*, int&, double*,
		int&, double&, double*, int&);

using namespace std;

const string PolynomialSurface::name = "Polynomial";

// ____________________________________________________________________________
// Creation, Destruction, Initialization 
// ____________________________________________________________________________


PolynomialSurface::PolynomialSurface(SurfData& sd, unsigned order, 
  unsigned responseIndex) : Surface(0),  
  xsize(sd.xSize()), order(order), digits(order) 
{
#ifdef __TESTING_MODE__
  constructCount++;
#endif
  dataItr = new SurfDataIterator(sd, responseIndex);
  resetTermCounter();
  build();
}  

PolynomialSurface::PolynomialSurface(AbstractSurfDataIterator* dataItr, 
  unsigned order) : Surface(dataItr), xsize(dataItr->xSize()), order(order),
  digits(order)
{
#ifdef __TESTING_MODE__
  constructCount++;
#endif
  resetTermCounter();
  build();
}
  
PolynomialSurface::PolynomialSurface(unsigned xsize, unsigned order, 
  std::vector<double> coefficients) : Surface(0), 
  xsize(xsize), order(order), coefficients(coefficients),
  digits(order)
{
#ifdef __TESTING_MODE__
  constructCount++;
#endif
  resetTermCounter();
  build();
}

PolynomialSurface::PolynomialSurface(const string filename) : Surface(0)
{
#ifdef __TESTING_MODE__
  constructCount++;
#endif
  read(filename);
}

PolynomialSurface::~PolynomialSurface()
{
#ifdef __TESTING_MODE__
  destructCount++;
#endif
}

// ____________________________________________________________________________
// Queries 
// ____________________________________________________________________________

const string PolynomialSurface::surfaceName() const
{
  return name;
}

unsigned PolynomialSurface::minPointsRequired(unsigned xsize, unsigned order)
{
  return nChooseR(xsize + order, order);
}

unsigned PolynomialSurface::minPointsRequired() const
{ 
  return minPointsRequired(xsize, order);
}

double PolynomialSurface::evaluate(const std::vector<double> & x)
{ 
  double sum = 0;
  resetTermCounter();
  for (unsigned i = 0; i < coefficients.size(); i++) {
    sum += coefficients[i] * computeTerm(x);
    nextTerm();
  }
  return sum;
}

//double PolynomialSurface::errorMetric(string metricName)
//{
//  if (metricName == "press") {
//    //cout << "Call Press from Polynomial" << endl;
//    return press();
//  } else if (metricName == "rsquared") {
//    //cout << "Call rSquared from Polynomial" << endl;
//    return rSquared();
//  } else {
//    return Surface::errorMetric(metricName);	
//  }
//  return 0;
//}

//double PolynomialSurface::press()
//{
//  cout << "Polynomial::press" << endl;
//  return 0;
//}

//double PolynomialSurface::rSquared(AbstractSurfDataIterator* iter)
//{
//  if(!iter) {
//    iter = this->dataItr;
//    if (!iter) {
//      cerr << "Cannot compute R^2 without data" << endl;
//      return 0.0;
//    }
//  }
//
//  int pts = iter->elementCount();
//  int lwork = pts * coefficients.size() * 2;
//
//  // allocate space for the matrices
//  double* a = new double[coefficients.size()*pts];
//  double* a2 = new double[coefficients.size()*pts];
//  double* b = new double[pts];
//  double* y = new double[pts];
//  double* yhat = new double[pts];
//  double* work = new double[lwork];
//
//  double ysum = 0.0;
//  // populate the A and B matrices in preparation for Ax=b
//  iter->toFront();
//  unsigned point = 0;
//  while(!iter->isEnd()) {
//    resetTermCounter();
//    while (termIndex < coefficients.size()) {
//          a[point+termIndex*pts] = computeTerm(iter->currentElement().X());
//          a2[point+termIndex*pts] = computeTerm(iter->currentElement().X());
//          //cout << "a[" << point+termIndex*pts << "] = " << a[point+termIndex*pts] << endl;
//        nextTerm();
//    }
//    b[point] = iter->currentElement().F(iter->responseIndex());
//    y[point] = b[point];
//    ysum += y[point];
//    iter->nextElement();
//    point++;
//  }
//  double ysumSquared = ysum * ysum;
//  int inc = 1;
//  double ydoty = ddot_(pts,y,inc,y,inc);
//  // values must be passed by reference to Fortran, so variables must be declared for info, nrhs, trans
//  int info;
//  int nrhs=1;
//  char trans = 'N';
//  int numCoeff = static_cast<int>(coefficients.size());
//  //cout << "A Matrix: " << endl;
//  //writeMatrix(a,pts,numCoeff,cout);
//  //cout << "B Vector: " << endl;
//  //long defaultPrecision = cout.precision();
//  //cout.precision(30);
//  //writeMatrix(b,pts,1,cout);
//  //cout.precision(defaultPrecision);
//  dgels_(trans,pts,numCoeff,nrhs,a,pts,b,pts,work,lwork,info);
//  //cout << "A Matrix after: " << endl;
//  //writeMatrix(a,pts,numCoeff,cout);
//  if (info < 0) {
//          cerr << "dgels_ returned with an error" << endl;
//  }
//  double noScale = 1.0;
//  double noExist = 0.0;
//  dgemv_(trans,pts,numCoeff,noScale,a2,pts,b,inc,noExist,yhat,inc);
//  //cout << "Yhat vector: " << endl;
//  //cout.precision(30);
//  //writeMatrix(yhat,pts,1,cout);
//  //cout.precision(defaultPrecision);
//  double yhatDoty = ddot_(pts,yhat,inc,y,inc);
//  double Syy = ydoty - ysumSquared;
//  double SSe = ydoty - yhatDoty;
//  double r2 = 1.0 - SSe/Syy;
//  
//  delete [] work;
//  delete [] b;
//  delete [] a;
//  delete [] a2;
//  delete [] yhat;
//  delete [] y;
//  cout << "ysumSquared: " << ysumSquared << " ydoty: " << ydoty << endl;
//  return r2;
//}

// ____________________________________________________________________________
// Commands 
// ____________________________________________________________________________

void PolynomialSurface::build()
{
  if (!acceptableData()) {
    cerr << "Cannot build surface.  Data are not acceptable" << endl;
  } else {
    coefficients.resize(minPointsRequired());
    int pts = dataItr->elementCount();
    int lwork = pts * coefficients.size() * 2;

    // allocate space for the matrices
    double* a = new double[coefficients.size()*pts];
    double* b = new double[pts];
    double* work = new double[lwork];

    // populate the A and B matrices in preparation for Ax=b
    dataItr->toFront();
    unsigned point = 0;
    while(!dataItr->isEnd()) {
      resetTermCounter();
      while (termIndex < coefficients.size()) {
        a[point+termIndex*pts] = computeTerm(dataItr->currentElement().X());
        nextTerm();
      }
      b[point] = dataItr->currentElement().F(dataItr->responseIndex());
      dataItr->nextElement();
      point++;
    }

    // values must be passed by reference to Fortran, so variables must be declared for info, nrhs, trans
    int info;
    int nrhs=1;
    char trans = 'N';
    int numCoeff = static_cast<int>(coefficients.size());
    writeMatrix("AMatrix",a,static_cast<unsigned>(pts),coefficients.size());
    writeMatrix("BVector",b,static_cast<unsigned>(pts),1);
    dgels_(trans,pts,numCoeff,nrhs,a,pts,b,pts,work,lwork,info);
    //cout << "A Matrix after: " << endl;
    //writeMatrix(a,pts,numCoeff,cout);
    if (info < 0) {
            cerr << "dgels_ returned with an error" << endl;
    }

    for(unsigned i=0;i<coefficients.size();i++) {
      coefficients[i] = b[i];
      double approx = floor(coefficients[i]+.5);
      if (abs(approx-coefficients[i]) < 1.0e-5) {
      	coefficients[i] = approx;
      }
    }

    //writeText(cout);
    delete [] work;
    delete [] b;
    delete [] a;
  }
}

/// Create a surface of the same type as 'this.'  This objects data should
/// be replaced with the dataItr passed in, but all other attributes should
/// be the same (e.g., a second-order polynomial should return another 
/// second-order polynomial.  Surfaces returned by this method can be used
/// to compute the PRESS statistic.
PolynomialSurface* PolynomialSurface::makeSimilarWithNewData
  (AbstractSurfDataIterator* dataItr)
{
  return new PolynomialSurface(dataItr, order);
}

// ____________________________________________________________________________
// Helper methods 
// ____________________________________________________________________________

unsigned PolynomialSurface::fact(unsigned x)
{
  if (x > 12) {
    cerr << x << "! exceeds the size of an integer. Don't do it." << endl;
    return 1;
  }

  unsigned result = 1;
  for (unsigned i = 1; i <= x; i++) {
    result *= i;
  }
  return result;
}

unsigned PolynomialSurface::nChooseR(unsigned n, unsigned r)
{
  //return fact(n) / (fact(n-r)*fact(r));
  unsigned num = 1;
  unsigned den = 1;
  unsigned lim = (n-r) < r ? n-r : r;
  for (unsigned i = 0; i < lim; i++) {
    num *= (n-i);
    den *= i+1;
  }
  return num/den;
}
	   
void PolynomialSurface::resetTermCounter()
{
  for (unsigned i = 0; i < digits.size(); i++) {
    digits[i] = 0;
  }
  termIndex = 0;
}

double PolynomialSurface::computeTerm(const std::vector<double>& x)
{
  double product = 1.0;
  for (unsigned i = 0; i < digits.size(); i++) {
    if (digits[i] != 0) {
      if (digits[i] - 1 >= x.size()) {
	cerr << "dimension mismatch in PolynomialSurface::computeTerm" << endl;
	return 0.0;
      }
      product *= x[digits[i]-1];
    }
  }
  return product;
}

void PolynomialSurface::nextTerm()
{
  unsigned curDig = 0;
  while (digits[curDig] == xsize && curDig < order) {
    curDig++;
  }
  // Don't go past last term
  if (curDig < order) {
    digits[curDig]++;
    while(curDig > 0) {
      curDig--;
      digits[curDig] = digits[curDig+1];
    }
  }
  termIndex++;
}

// ____________________________________________________________________________
// I/O 
// ____________________________________________________________________________

/// Save the surface to a file in binary format
void PolynomialSurface::writeBinary(ostream& os)
{
  os.write(reinterpret_cast<char*>(&xsize),sizeof(xsize));
  os.write(reinterpret_cast<char*>(&order),sizeof(order));
  for (unsigned i = 0; i < coefficients.size(); i++) {
    os.write(reinterpret_cast<char*>(&coefficients[i]),sizeof(coefficients[i]));
  }
  // figure out what to do with data
}

void PolynomialSurface::writeText(ostream& os) 
{
  resetTermCounter();
  os << xsize << " dimensions" << endl;
  os << order << " order" << endl;
  os.setf(ios::left);
  for (unsigned i = 0; i < coefficients.size(); i++) {
    os << setprecision(17) << setw(26) << coefficients[i];
    printTermLabel(os);
    nextTerm();
    if (i+1 < coefficients.size()) {
      os << " +";
    }
    os << endl;
  }
  os.unsetf(ios::left);
}

/// Load the surface from a file
void PolynomialSurface::readBinary(istream& is)
{
  is.read(reinterpret_cast<char*>(&xsize),sizeof(xsize));
  is.read(reinterpret_cast<char*>(&order),sizeof(order));
  coefficients.resize(minPointsRequired());
  for (unsigned i = 0; i < coefficients.size(); i++) {
    is.read(reinterpret_cast<char*>(&coefficients[i]),sizeof(coefficients[i]));
  }
  digits.resize(order);
  resetTermCounter();
  valid = true;
  //originalData = false;
}

void PolynomialSurface::readText(istream& is)
{
  // read in the number of dimensions
  string sline;
  getline(is,sline);   
  istringstream streamline(sline);
  streamline.str(sline);
  streamline >> xsize; 		 
  // read in the order of the model
  getline(is,sline);   
  streamline.str(sline);
  streamline >> order;
  // determine the number of terms that should be in the file
  coefficients.resize(minPointsRequired());
  for (unsigned i = 0; i < coefficients.size(); i++) {
    getline(is,sline);   
    streamline.str(sline);
    // read each coefficient; ignore the label, if any, to the right 
    streamline >> coefficients[i];		
  }
  digits.resize(order);
  resetTermCounter();
  valid = true;
  //originalData = false;
}

void PolynomialSurface::printTermLabel(ostream& os)
{
  ostringstream labelStream;
  bool needStar = false;
  unsigned i = 0;
  while (i < order) {
	if (digits[i] != 0) {
          unsigned var = digits[i];
	    unsigned count = 1;
	    while (i+1 < order && digits[i+1] == var) {
		i++;
		count++;
	    }
	    if (needStar) {
	        labelStream << "*";
	    } else {
		needStar = true;
	    }
	    labelStream << "x" << var;
	    if (count > 1) {
	        labelStream << "^" << count;
	    }	
	}
	i++;
  }
  string label = labelStream.str();
  os << label;
}

  
// ____________________________________________________________________________
// Testing 
// ____________________________________________________________________________

#ifdef __TESTING_MODE__
  int PolynomialSurface::constructCount = 0;
  int PolynomialSurface::copyCount = 0;
  int PolynomialSurface::destructCount = 0;
#endif

