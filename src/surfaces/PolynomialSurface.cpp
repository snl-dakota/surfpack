#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include "PolynomialSurface.h"
#include "SurfData.h"

#include <cmath>
#include <iostream>
extern "C" void dgels_(char&, int&, int&, int&, double*,
                       int&, double*, int&, double*, int&, int&);

extern "C" double ddot_(int&, double*, int&, double*, int&);

extern "C" void dgemv_(char&, int&, int&, double&, double*, int&, double*,
		int&, double&, double*, int&);

using namespace std;

//PolynomialSurface::PolynomialSurface(const std::vector<double> & coeff)
//   : Surface("PolynomialSurface",coeff) { needsRebuild = false; }

PolynomialSurface::PolynomialSurface(string filename) : Surface(0, -1), k(0), n(0), numCoeff(0)
{
    loadBinary(filename);
}

PolynomialSurface::PolynomialSurface(int numvars, int order, std::vector<double> coefficients)
   : Surface(0, -1), k(order), n(numvars), numCoeff(coefficients.size())
{
    digits.resize(k);
    resetTermCounter();
    this->coefficients = coefficients;
}

PolynomialSurface::PolynomialSurface(istream& is, int responseIndex)
   : Surface(0, responseIndex)
{
    const int MAX_CHAR = 200;
    char line[MAX_CHAR];
    is.getline(line, MAX_CHAR);    // throw away first line; header info only
    is.getline(line, MAX_CHAR);   
    string sline(line);
    istringstream streamline(sline);
    streamline >> k; 		// kth order polynomial
    is.getline(line, MAX_CHAR);   
    sline = line;
    streamline.str(sline);
    streamline >> n;		// n variables
    is.getline(line, MAX_CHAR);   
    sline = line;
    streamline.str(sline);
    streamline >> numCoeff;		// number of coefficients
    coefficients.resize(numCoeff);
    for (int i = 0; i < numCoeff; i++) {
        is.getline(line, MAX_CHAR);   
        sline = line;
        streamline.str(sline);
        streamline >> coefficients[i];		// read each coefficient; ignore the label, if any, to the right 
    }
    digits.resize(k);
    resetTermCounter();
    needsRebuild = false;
}

PolynomialSurface::PolynomialSurface(SurfData * data, int order, int responseIndex)
   : Surface(data, responseIndex), k(order) 
{ 
    n = data->getDimension();
    digits.resize(k);
    resetTermCounter();
    numCoeff = PolynomialSurface::getMinPointCount(k,n);
    responseIndex = 0;

}

/// Save the surface to a file in binary format
void PolynomialSurface::saveBinary(std::string filename)
{
    ofstream outfile(filename.c_str(),ios::out | ios::binary);
    if (!outfile) {
	    cerr << "Could not open " << filename << " for writing." << endl;
	    return;
    }
	
    int tempID = Surface::polynomialSurfaceID;
    outfile.write((char*)(&tempID),sizeof(int));
    outfile.write((char*)(&k),sizeof(int));
    outfile.write((char*)(&n),sizeof(int));
    outfile.write((char*)(&numCoeff),sizeof(int));
    for (int i = 0; i < numCoeff; i++) {
	outfile.write((char*)(&coefficients[i]),sizeof(double));
    }
    surfData->write(outfile, true); // binary write
    outfile.close();
}

/// Load the surface from a file
void PolynomialSurface::loadBinary(std::string filename)
{
    ifstream infile(filename.c_str(),ios::in | ios::binary);
    if (!infile) {
	    cerr << "Could not open " << filename << " for reading." << endl;
	    return;
    }
	
    int tempID; // = Surface::polynomialSurfaceID;
    infile.read((char*)(&tempID),sizeof(int));
    if (tempID != Surface::polynomialSurfaceID) {
	    cerr << "ID in surface file does not specify Polynomial" << endl;
	    return;
    }
    infile.read((char*)(&k),sizeof(int));
    infile.read((char*)(&n),sizeof(int));
    infile.read((char*)(&numCoeff),sizeof(int));
    coefficients.resize(numCoeff);
    for (int i = 0; i < numCoeff; i++) {
	infile.read((char*)(&coefficients[i]),sizeof(double));
    }
    delete surfData;
    surfData = new SurfData;
    surfData->read(infile,true); // binary read
    //surfData->write(cout); // text write
    infile.close();
    digits.resize(k);
    resetTermCounter();
    needsRebuild = false;
}

void PolynomialSurface::save(std::string filename)
{
    ofstream outfile(filename.c_str(),ios::out);
    if (!outfile) {
	    cerr << "Could not open " << filename << " for writing." << endl;
	    return;
    }

    outfile << "Polynomial Surface Fit" << endl;
    outfile << k << "   order" << endl
	    << n << "   variables" << endl
	    << numCoeff << "   number of coefficients" << endl;
    write(outfile);
    outfile.close();
}

//PolynomialSurface::PolynomialSurface(const PolynomialSurface & surf)
//   : Surface(surf) { needsRebuild = false; }

PolynomialSurface::~PolynomialSurface() {}



int PolynomialSurface::fact(int x)
{
   if (x > 12) {
	   cerr << x << "! exceeds the size of an integer. Don't do it." << endl;
	   return 1;
   }
   int result = 1;
   for (int i = 1; i <= x; i++) {
       result *= i;
   }
   return result;
}

int PolynomialSurface::nChooseR(int n, int r)
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
	current = 0;
}

//void print()
//{
//    for (vector<int>::reverse_iterator it = digits.rbegin(); it != digits.rend(); ++it) {
//        cout << *it;
//    }
//    cout << "  "; // endl;
//}

void PolynomialSurface::printTermLabel(ostream& os)
{
    ostringstream labelStream;
    bool needStar = false;
    int i = 0;
    while (i < k) {
	if (digits[i] != 0) {
            int var = digits[i];
	    int count = 1;
	    while (i+1 < k && digits[i+1] == var) {
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

double PolynomialSurface::computeTerm(const std::vector<double>& x)
{
    double product = 1.0;
    for (unsigned i = 0; i < digits.size(); i++) {
        if (digits[i] != 0) {
	    if (static_cast<unsigned>(digits[i]) - 1 >= x.size()) {
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
    int curDig = 0;
    while (digits[curDig] >= n && curDig < k) {
        curDig++;
    }
    if (curDig < k) {
	digits[curDig]++;
	while(curDig > 0) {
	    curDig--;
	    digits[curDig] = digits[curDig+1];
	}
    }
    current++;
}

//int main() 
//{
//	numCoeff = nChooseR(n+k,k);
//	reset();
//	do {
//	    print();
//	    printLabel();
//	    cout << endl;
//	    next();
//	} while (current < numCoeff);
//	cout << "numCoefficients: " << numCoeff << endl;
//
//	return 0;
//}
//


int PolynomialSurface::getMinPointCount(int order, int numvars)
{
    return nChooseR(order+numvars, order);
}

int PolynomialSurface::getMinPointCount(int val) const
{ 
    int result = PolynomialSurface::getMinPointCount(val,n);
    return result;
}

int PolynomialSurface::getDimension() const
{
    return n; // number of dimensions set in constructor
}
   
//double PolynomialSurface::getError(int code) const
//{ return 0; }

double PolynomialSurface::calculate(const std::vector<double> & x)
   throw(SurfException)
{ 
    double sum = 0;
    resetTermCounter();
    for (int i = 0; i < numCoeff; i++) {
        sum += coefficients[i] * computeTerm(x);
	nextTerm();
    }
    return sum;
}
    
void PolynomialSurface::calculateInternals(AbstractSurfDataIterator* iter)
{
   if(surfData==0) return;

   // m = pts
   // n = coeffCnt

   //int dim = surfData->getDimension();
   int pts = iter->getElementCount();
   int lwork = pts * numCoeff * 2;

   // allocate space for the matrices
   double* a = new double[numCoeff*pts];
   double* b = new double[pts];
   double* work = new double[lwork];

   // populate the A and B matrices in preparation for Ax=b
   iter->toFront();
   int point = 0;
   while(!iter->isEnd()) {
       resetTermCounter();
       while (current < numCoeff) {
	   a[point+current*pts] = computeTerm(iter->currentElement()->getX());
	   //cout << "a[" << point+current*pts << "] = " << a[point+current*pts] << endl;
           nextTerm();
       }
       b[point] = iter->currentElement()->getF(responseIndex);
       iter->nextElement();
       point++;
   }

   // values must be passed by reference to Fortran, so variables must be declared for info, nrhs, trans
   int info;
   int nrhs=1;
   char trans = 'N';
   //cout << "A Matrix: " << endl;
   //writeMatrix(a,pts,numCoeff,cout);
   //cout << "B Vector: " << endl;
   //writeMatrix(b,pts,1,cout);
   dgels_(trans,pts,numCoeff,nrhs,a,pts,b,pts,work,lwork,info);
   //cout << "A Matrix after: " << endl;
   //writeMatrix(a,pts,numCoeff,cout);
   if (info < 0) {
	   cerr << "dgels_ returned with an error" << endl;
   }

   coefficients.clear();
   coefficients.resize(numCoeff);

   for(int i=0;i<numCoeff;i++) {
       coefficients[i] = b[i];
       double approx = floor(coefficients[i]+.5);
       if (abs(approx-coefficients[i]) < 1.0e-5) {
       	coefficients[i] = approx;
       }
   }

   delete [] work;
   delete [] b;
   delete [] a;
}

double PolynomialSurface::rSquared(AbstractSurfDataIterator* iter)
{
   if(!surfData) {
	  cerr << "Cannot compute R^2 without data" << endl;
	  return 0.0;
   }

   bool deleteIter = false;
   if (!iter) {
	   iter = new SurfDataIterator(surfData);
	   deleteIter = true;
   }

   
   int pts = iter->getElementCount();
   int lwork = pts * numCoeff * 2;

   // allocate space for the matrices
   double* a = new double[numCoeff*pts];
   double* a2 = new double[numCoeff*pts];
   double* b = new double[pts];
   double* y = new double[pts];
   double* yhat = new double[pts];
   double* work = new double[lwork];

   double ysum = 0.0;
   // populate the A and B matrices in preparation for Ax=b
   iter->toFront();
   int point = 0;
   while(!iter->isEnd()) {
       resetTermCounter();
       while (current < numCoeff) {
	   a[point+current*pts] = computeTerm(iter->currentElement()->getX());
	   a2[point+current*pts] = computeTerm(iter->currentElement()->getX());
	   //cout << "a[" << point+current*pts << "] = " << a[point+current*pts] << endl;
           nextTerm();
       }
       b[point] = iter->currentElement()->getF(responseIndex);
       y[point] = b[point];
       ysum += y[point];
       iter->nextElement();
       point++;
   }
   double ysumSquared = ysum * ysum;
   int inc = 1;
   double ydoty = ddot_(pts,y,inc,y,inc);
   // values must be passed by reference to Fortran, so variables must be declared for info, nrhs, trans
   int info;
   int nrhs=1;
   char trans = 'N';
   //cout << "A Matrix: " << endl;
   //writeMatrix(a,pts,numCoeff,cout);
   //cout << "B Vector: " << endl;
   //long defaultPrecision = cout.precision();
   //cout.precision(30);
   //writeMatrix(b,pts,1,cout);
   //cout.precision(defaultPrecision);
   dgels_(trans,pts,numCoeff,nrhs,a,pts,b,pts,work,lwork,info);
   //cout << "A Matrix after: " << endl;
   //writeMatrix(a,pts,numCoeff,cout);
   if (info < 0) {
	   cerr << "dgels_ returned with an error" << endl;
   }
   double noScale = 1.0;
   double noExist = 0.0;
   dgemv_(trans,pts,numCoeff,noScale,a2,pts,b,inc,noExist,yhat,inc);
   //cout << "Yhat vector: " << endl;
   //cout.precision(30);
   //writeMatrix(yhat,pts,1,cout);
   //cout.precision(defaultPrecision);
   double yhatDoty = ddot_(pts,yhat,inc,y,inc);
   double Syy = ydoty - ysumSquared;
   double SSe = ydoty - yhatDoty;
   double r2 = 1.0 - SSe/Syy;
   
   delete [] work;
   delete [] b;
   delete [] a;
   delete [] a2;
   delete [] yhat;
   delete [] y;
   cout << "ysumSquared: " << ysumSquared << " ydoty: " << ydoty << endl;
   if (deleteIter) {
	   delete iter;
   }
   return r2;
}


//ostream & PolynomialSurface::writeClean(ostream & os) throw(SurfException)
//{
//    resetTermCounter();
//    os.setf(ios::left);
//    for (int i = 0; i < numCoeff; i++) {
//        os << setprecision(25) << setw(35) << coefficients[i];
//	printTermLabel(os);
//	nextTerm();
//	if (i+1 < numCoeff) {
//		os << " +";
//	}
//	os << endl;
//    }
//    os.unsetf(ios::left);
//    return os;
//}

ostream& PolynomialSurface::write(ostream& os) 
{
    resetTermCounter();
    os.setf(ios::left);
    for (int i = 0; i < numCoeff; i++) {
        os << setprecision(25) << setw(35) << coefficients[i];
	printTermLabel(os);
	nextTerm();
	if (i+1 < numCoeff) {
		os << " +";
	}
	os << endl;
    }
    os.unsetf(ios::left);
    return os;
    //return writeClean(os);
}

string PolynomialSurface::getType() const
{
    return "Polynomial";
}
   
ostream& PolynomialSurface::writeMatrix(double* mat, int rows, int columns, ostream& os)
{
    for (int r = 0; r < rows; r++) {
	for (int c = 0; c < columns; c++) {
	    os << setw(15) << mat[r + c*rows];
	}
	os << endl;
    }
    return os;
}

double PolynomialSurface::errorMetric(string metricName)
{
	if (metricName == "press") {
		cout << "Call Press from Polynomial" << endl;
		return computePressStatistic();
	} else if (metricName == "rsquared") {
		cout << "Call rSquared from Polynomial" << endl;
		return rSquared();
	} else {
		return Surface::errorMetric(metricName);	
	}
	return 0;

}

double PolynomialSurface::press()
{
	cout << "Polynomial::press" << endl;
	return 0;

}

double PolynomialSurface::rsquared()
{
	cout << "Polynomial::rsquared" << endl;
	return 0;

}
