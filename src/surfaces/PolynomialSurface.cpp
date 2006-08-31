/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#include "surfpack_config.h"
#include "surfpack.h"
#include "SurfData.h"
#include "PolynomialSurface.h"

using std::endl;
using std::ios;
using std::istream;
using std::istringstream;
using std::ostream;
using std::ostringstream;
using std::range_error;
using std::setprecision;
using std::setw;
using std::string;
using std::vector;

const string PolynomialSurface::name = "polynomial";

// ____________________________________________________________________________
// Creation, Destruction, Initialization 
// ____________________________________________________________________________


PolynomialSurface::PolynomialSurface(SurfData* sd, unsigned order)
  : Surface(sd), order(order), digits(order) 
{
#ifdef __TESTING_MODE__
  constructCount++;
#endif
  resetTermCounter();
}  

PolynomialSurface::PolynomialSurface(unsigned x_size, unsigned order, 
  vector<double> coefficients) : Surface(0), order(order), 
  coefficients(coefficients), digits(order)
{
#ifdef __TESTING_MODE__
  constructCount++;
#endif
  this->xsize = x_size;
  builtOK = true;
  //dataAdded = false;
  //dataModified = false;
  resetTermCounter();
}

PolynomialSurface::PolynomialSurface(const string filename) : Surface(0)
{
#ifdef __TESTING_MODE__
  constructCount++;
#endif
  read(filename);
}

PolynomialSurface::PolynomialSurface(const PolynomialSurface& other) 
  : Surface(other), order(other.order), coefficients(other.coefficients),
  digits(other.order), termIndex(other.termIndex) 
{
#ifdef __TESTING_MODE__
  constructCount++;
#endif
  // Must assign after base constructor has been called, because if the
  // surface has been built, the setData method called in the base
  // copy constructor will set dataModified to true, perhaps unnecessarily
  dataModified = other.dataModified;
}

/// Create a surface of the same type as 'this.'  This objects data should
/// be replaced with the dataItr passed in, but all other attributes should
/// be the same (e.g., a second-order polynomial should return another 
/// second-order polynomial.  Surfaces returned by this method can be used
/// to compute the PRESS statistic.
PolynomialSurface* 
  PolynomialSurface::makeSimilarWithNewData(SurfData* surfData)
{
  PolynomialSurface* newPS = new PolynomialSurface(*this);
  newPS->setData(surfData);
  newPS->createModel();
  return newPS;
}

PolynomialSurface::~PolynomialSurface()
{
#ifdef __TESTING_MODE__
  destructCount++;
#endif
}

void PolynomialSurface::config(const Arg& arg)
{
  string argname = arg.name;
  if (argname == "order") {
    setOrder(arg.getRVal()->getInteger());
  } else {
    Surface::config(arg);
  }
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
  //return nChooseR(xsize + order, order);
  vector<double> coeff;
  PolynomialSurface tempps(xsize,order,coeff);
  tempps.eqConRHS.clear(); // Make sure it's empty
  return tempps.minPointsRequired();
}

unsigned PolynomialSurface::minPointsRequired() const
{ 
  if (xsize == 0) {
    throw string("Cannot create surface with zero dimensionality");
  } else if (order == 0) {
    return 1;
  } else {
    //unsigned result = minPointsRequired(xsize, order);
    //return result; 
    resetTermCounter();
    while (!lastTerm) {
      nextTerm();
    }
    // Add one because the terms are numbered 0 to n-1
    return termIndex + 1 - eqConRHS.size();
  }
}

double PolynomialSurface::evaluate(const vector<double> & x)
{ 
  // Add sanity check on x
  double sum = 0;
  resetTermCounter();
  for (unsigned i = 0; i < coefficients.size(); i++) {
    sum += coefficients[i] * computeTerm(x);
    nextTerm();
  }
  return sum;
}

void PolynomialSurface::gradient(const vector<double> & x, vector<double>& gradient_vector)
{ 
  // Add sanity check on x
  assert(x.size() == xsize);
  gradient_vector = vector<double>(xsize,0.0);
  vector<unsigned> factorCounts;
  // Entry difVar tells how many times to differentiate with respect to the ith var
  // For example d^2f/dx1^2 would have a 2 for the x1 entry, d^2f/dx1,x2 would
  // have 1s for both the x1 and x2 entries.
  vector<unsigned> differentiationCounts;
  // Differentiate with respect to each variable
  for (unsigned difVar = 0; difVar < xsize; difVar++) {
    resetTermCounter();
    differentiationCounts = vector<unsigned>(xsize,0);
    differentiationCounts[difVar] = 1;
    for (unsigned term = 0; term < coefficients.size(); term++) {
      accumulateLikeFactors(factorCounts);
      gradient_vector[difVar] += coefficients[term] * computeDerivTerm(x, factorCounts, differentiationCounts);
      nextTerm();
    }
  }
}

void PolynomialSurface::hessian(const vector<double> & x, SurfpackMatrix<double>& hessian)
{ 
  // Add sanity check on x
  assert(x.size() == xsize);
  hessian.reshape(xsize,xsize);
  vector<unsigned> factorCounts;
  // Entry difVar tells how many times to differentiate with respect to the ith var
  // For example d^2f/dx1^2 would have a 2 for the x1 entry, d^2f/dx1,x2 would
  // have 1s for both the x1 and x2 entries.
  vector<unsigned> differentiationCounts;
  // Differentiate with respect to each variable
  for (unsigned difVar1 = 0; difVar1 < xsize; difVar1++) {
    for (unsigned difVar2 = difVar1; difVar2 < xsize; difVar2++) {
      hessian[difVar1][difVar2] = hessian[difVar2][difVar1] = 0.0;
      resetTermCounter();
      differentiationCounts = vector<unsigned>(xsize,0);
      differentiationCounts[difVar1]++;
      differentiationCounts[difVar2]++;
      for (unsigned term = 0; term < coefficients.size(); term++) {
        accumulateLikeFactors(factorCounts);
        double addition = coefficients[term] * computeDerivTerm(x, factorCounts, differentiationCounts);
        hessian[difVar1][difVar2] += addition;
        if (difVar1 != difVar2) hessian[difVar2][difVar1] += addition;
        nextTerm();
      }
    }
  }
}

void PolynomialSurface::setEqualityConstraints(unsigned asv,const SurfPoint& sp,  double valuePtr, vector<double>* gradientPtr, SurfpackMatrix<double>* hessianPtr)
{
  coefficients.resize(minPointsRequired(xsize,order));
  unsigned numConstraints = 0;
  if (asv & 1) numConstraints += 1; // value at a particular point
  if (asv & 2) numConstraints += xsize; // gradient at a point
  if (asv & 4) numConstraints += (xsize*xsize+xsize)/2; // hessian at a point
  eqConRHS.resize( numConstraints );
  // Must compute number of terms first
  SurfpackMatrix<double> temp(eqConRHS.size(),coefficients.size(),true);
  eqConLHS = temp;
  //eqConLHS.reshape(eqConRHS.size(),coefficients.size());
  // Marks the index of the next constraint to be added (necessary since
  // indices of e.g. the gradient constraints will be different depending on
  // whether or not the value constraint is used
  unsigned index = 0;
  // If requested, add the equality constraint for the point value
  if (asv & 1) {
    resetTermCounter();
    while (!lastTerm) {
      eqConLHS[index][termIndex] = computeTerm(sp.X());
      nextTerm();
    }
    eqConRHS[index] = valuePtr;
    ++index;
  }

  // If requested, add the equality constraints for the gradient
  if (asv & 2) {
    const vector<double>& gradient = *gradientPtr;
    assert(gradient.size() == xsize);
    vector<unsigned> factorCounts;
    vector<unsigned> differentiationCounts;
    for (unsigned difVar = 0; difVar < xsize; difVar++ ) {
      differentiationCounts = vector<unsigned>(xsize,0);
      differentiationCounts[difVar] = 1;
      resetTermCounter(); 
      while (!lastTerm) {
        accumulateLikeFactors(factorCounts);
        eqConLHS[index][termIndex] = 
          computeDerivTerm(sp.X(), factorCounts, differentiationCounts);
        nextTerm();
      }
      eqConRHS[index] = gradient[difVar];
      ++index;
    }
  }

  // If requested, add the equality constraints for the hessian
  if (asv & 4) {
    SurfpackMatrix<double>& hessian = *hessianPtr;
    assert(hessian.getNCols() == xsize);
    assert(hessian.getNRows() == xsize);
    vector<unsigned> factorCounts;
    vector<unsigned> differentiationCounts;
    for (unsigned difVar1 = 0; difVar1 < xsize; difVar1++ ) {
      for (unsigned difVar2 = difVar1; difVar2 < xsize; difVar2++ ) {
        differentiationCounts = vector<unsigned>(xsize,0);
        differentiationCounts[difVar1]++;
        differentiationCounts[difVar2]++;
        resetTermCounter(); 
        while (!lastTerm) {
          accumulateLikeFactors(factorCounts);
          eqConLHS[index][termIndex] = 
            computeDerivTerm(sp.X(), factorCounts, differentiationCounts);
          nextTerm();
        }
        eqConRHS[index] = hessian[difVar1][difVar2];
        ++index;
      } // difVar2
    } // difVar1
  } // if hessian needed
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
//  dgels(trans,pts,numCoeff,nrhs,a,pts,b,pts,work,lwork,info);
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

PolynomialSurface& PolynomialSurface::operator=(const PolynomialSurface& other)
{
  ///\todo check for assignment to self
  order = other.order;
  coefficients = other.coefficients;
  digits = other.digits;
  eqConLHS = other.eqConLHS;
  eqConRHS = other.eqConRHS;
  termIndex = other.termIndex;
  lastTerm = other.lastTerm;
  return *this;
}

void PolynomialSurface::build(SurfData& data)
{
  coefficients.resize(static_cast<int>(minPointsRequired(xsize,order)));
  SurfpackMatrix<double> A(data.size(),coefficients.size(),true);
  vector<double> b(data.size());
  // populate the A and B matrices in preparation for Ax=b
  for(unsigned i = 0; i < data.size(); i++) {
    resetTermCounter();
    while (!lastTerm) {
      A[i][termIndex] = computeTerm(data[i].X());
      nextTerm();
    }
    b[i] = data.getResponse(i);
  }
  // Solve the system of equations
  if (eqConRHS.empty()) {
    surfpack::linearSystemLeastSquares(A,coefficients,b);
  } else {
    //cout << "leastSquares with " << eqConRHS.size() << " equality constraints" << endl;
    surfpack::leastSquaresWithEqualityConstraints
      (A,coefficients,b,eqConLHS,eqConRHS);
  }
  // The mean of the predictions at the sample sites 
  // should be the same as the mean of training data responses
  double sum_fi = 0.0;
  double sum_fhati = 0.0;
  for (unsigned i = 0; i < data.size(); i++) {
    sum_fi += data.getResponse(i);
    sum_fhati += evaluate(data[i].X());
  }
  ///\todo Print out a warning if sum_fi and sum_fhati are too different
}

/// Set the degree of the polynomial fit (e.g., linear=1, quadratic=2, etc.)
void PolynomialSurface::setOrder(unsigned order_in)
{
    this->order = order_in;
    digits.resize(this->order);
}
// ____________________________________________________________________________
// Helper methods 
// ____________________________________________________________________________

//unsigned PolynomialSurface::fact(unsigned x)
//{
//  if (x > 12) {
//    cerr << x << "! exceeds the size of an integer. Don't do it." << endl;
//    return 1;
//  }
//
//  unsigned result = 1;
//  for (unsigned i = 1; i <= x; i++) {
//    result *= i;
//  }
//  return result;
//}

//unsigned PolynomialSurface::nChooseR(unsigned n, unsigned r)
//{
//  //return fact(n) / (fact(n-r)*fact(r));
//  unsigned num = 1;
//  unsigned den = 1;
//  unsigned lim = (n-r) < r ? n-r : r;
//  for (unsigned i = 0; i < lim; i++) {
//    num *= (n-i);
//    den *= i+1;
//  }
//  return num/den;
//}
	   
void PolynomialSurface::resetTermCounter() const
{
  for (unsigned i = 0; i < digits.size(); i++) {
    digits[i] = 0;
  }
  termIndex = 0;
  lastTerm = false;
}

double PolynomialSurface::computeTerm(const vector<double>& x) const
{
  double product = 1.0;
  for (unsigned i = 0; i < digits.size(); i++) {
    if (digits[i] != 0) {
      if (digits[i] - 1 >= x.size()) {
        ostringstream msg;
        msg << "Error in PolynomialSurface::computeTerm.  "
	    << "Variable " << digits[i] << " requested, but "
 	    << "point has only " << x.size() << " dimension(s)."
	    << endl;
        throw range_error(msg.str());
      }
      product *= x[digits[i]-1];
    }
  }
  return product;
}

double PolynomialSurface::computeDerivTerm(const vector<double>& x,
  const vector<unsigned>& factorCounts, 
  const vector<unsigned>& differentiationCounts) const
{
  double product = 1.0;
  for (unsigned var = 0; var < factorCounts.size(); var++) {
    if (factorCounts[var] < differentiationCounts[var])  {
      return 0.0;
    } else {
      for (unsigned i = 0; i < differentiationCounts[var]; i++) {
        product *= static_cast<double>(factorCounts[var]-i);
      }
      for (unsigned i = 0; i < factorCounts[var] - differentiationCounts[var];
		i++) {
        product *= x[var];
      }
    }
  } 
  return product;
}

void PolynomialSurface::accumulateLikeFactors(vector<unsigned>& factorCounts)
{
  factorCounts = vector<unsigned>(xsize,0);
  for (unsigned difVar = 0; difVar < digits.size(); difVar++) {
    if (digits[difVar]) factorCounts[digits[difVar]-1]++;
  }
}

void PolynomialSurface::nextTerm() const
{
  if (!lastTerm) {
    unsigned curDig = 0;
    while (curDig < order && digits[curDig] == xsize) {
      curDig++;
    }
    // Don't go past last term
    if (curDig < order) {
      digits[curDig]++;
      while(curDig > 0) {
        curDig--;
        digits[curDig] = digits[curDig+1];
      }
      if (termIndex == INT_MAX) {
        //for (unsigned j = 0; j < digits.size(); j++) {
        //  cout << digits[j] << endl;
        //}
        throw range_error(
          "Integer overflow: number of terms exceeds maximum integer");
      } else {
        termIndex++;
      }
    } else {
      lastTerm = true;
    }
  }
  //printTermLabel(cout);
  //printTermComponents(cout);
  //cout << endl;
 
}

// ____________________________________________________________________________
// I/O 
// ____________________________________________________________________________

/// Save the surface to a file in binary format
void PolynomialSurface::writeBinary(ostream& os)
{
  prepareData();
  os.write(reinterpret_cast<char*>(&xsize),sizeof(xsize));
  os.write(reinterpret_cast<char*>(&order),sizeof(order));
  for (unsigned difVar = 0; difVar < coefficients.size(); difVar++) {
    os.write(reinterpret_cast<char*>(&coefficients[difVar]),sizeof(coefficients[difVar]));
  }
  // figure out what to do with data
}

void PolynomialSurface::writeText(ostream& os) 
{
  prepareData();
  resetTermCounter();
  os << xsize << " dimension(s)" << endl;
  os << order << " order" << endl;
  os.setf(ios::left);
  for (unsigned difVar = 0; difVar < coefficients.size(); difVar++) {
    os << setprecision(17) << setw(26) << coefficients[difVar];
    printTermLabel(os);
    nextTerm();
    if (difVar+1 < coefficients.size()) {
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
  digits.resize(order);
  coefficients.resize(minPointsRequired(xsize,order));
  for (unsigned difVar = 0; difVar < coefficients.size(); difVar++) {
    is.read(reinterpret_cast<char*>(&coefficients[difVar]),sizeof(coefficients[difVar]));
  }
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
  digits.resize(order);
  // determine the number of terms that should be in the file
  coefficients.resize(minPointsRequired(xsize,order));
  for (unsigned difVar = 0; difVar < coefficients.size(); difVar++) {
    getline(is,sline);   
    streamline.str(sline);
    // read each coefficient; ignore the label, if any, to the right 
    streamline >> coefficients[difVar];		
  }
}

void PolynomialSurface::printTermLabel(ostream& os)
{
  ostringstream labelStream;
  bool needStar = false;
  unsigned difVar = 0;
  while (difVar < order) {
	if (digits[difVar] != 0) {
          unsigned var = digits[difVar];
	    unsigned count = 1;
	    while (difVar+1 < order && digits[difVar+1] == var) {
		difVar++;
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
	difVar++;
  }
  string label = labelStream.str();
  os << label;
}

void PolynomialSurface::printTermComponents(ostream& os)
{
  os << " ";
  for (unsigned difVar = 0; difVar < digits.size(); difVar++) {
    os << "[";
    if (digits[difVar] != 0) {
	os << digits[difVar];
    }
    os << "]";
  }
  os << " ";
}
  
// ____________________________________________________________________________
// Testing 
// ____________________________________________________________________________

#ifdef __TESTING_MODE__
  int PolynomialSurface::constructCount = 0;
  int PolynomialSurface::copyCount = 0;
  int PolynomialSurface::destructCount = 0;
#endif

