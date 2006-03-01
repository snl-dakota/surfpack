/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */
#include "surfpack_config.h"
#include "surfpack_system_headers.h"
#include "PolynomialSurfaceTest.h"
#include "SurfData.h"
#include "Surface.h"
#include "unittests.h"
#include "PolynomialSurface.h"

using namespace std;

CPPUNIT_TEST_SUITE_REGISTRATION( PolynomialSurfaceTest );

void PolynomialSurfaceTest::setUp()
{
  initialize();
}

void PolynomialSurfaceTest::tearDown()
{
}

//void PolynomialSurfaceTest::constructor()
//{
//  unsigned dim = 2500;
//  vector<SurfPoint> surfpoints;
//  surfpoints.reserve(dim);
//  cout << "Creating points" << endl;
//  for (unsigned pti = 0; pti < dim + 1; pti++) {
//    vector<double> pt(dim);
//    double response = 0.0;
//    for (unsigned dimi = 0; dimi < dim; dimi++) {
//      pt[dimi] = pti + dimi;
//      response += pt[dimi];
//    }
//    surfpoints.push_back(SurfPoint(pt, response));
//  }
//  cout << "Creating SurfData object" << endl;
//  SurfData sd(surfpoints); 
//  cout << "Creating PolynomialSurface object" << endl;
//  PolynomialSurface ps(&sd, 1);
//  ps.createModel();
//  cout << "Writing to file" << endl;
//  ps.write("manydim.txt");
//  CPPUNIT_ASSERT(1 == 1);
//}

//void PolynomialSurfaceTest::constructor()
//{
// 
//  unsigned dim = 5000;
//  vector<SurfPoint> surfpoints;
//  surfpoints.reserve(dim);
//  cout << "Creating points" << endl;
//  for (unsigned pti = 0; pti < dim + 1; pti++) {
//    vector<double> pt(1);
//    pt[0] = pti;
//    double response = pti;
//    surfpoints.push_back(SurfPoint(pt, response));
//  }
//  cout << "Creating SurfData object" << endl;
//  SurfData sd(surfpoints); 
//  cout << "Creating PolynomialSurface object" << endl;
//  PolynomialSurface ps(&sd, dim);
//  ps.createModel();
//  cout << "Writing to file" << endl;
//  ps.write("manydim.txt");
//  CPPUNIT_ASSERT(1 == 1);
//
//} 

//void PolynomialSurfaceTest::constructor()
//{
//  SurfData sd(fullPath("oneDimQuadratic.txt"));
//  PolynomialSurface ps(&sd, 2);
//  ps.createModel();
//  ps.write("oneDQpoly2.txt");
//  CPPUNIT_ASSERT(ps.coefficients.size() == 3);
//  CPPUNIT_ASSERT(matches(ps.coefficients[0],0.0));
//  CPPUNIT_ASSERT(matches(ps.coefficients[1],0.0));
//  CPPUNIT_ASSERT(matches(ps.coefficients[2],1.0));
//}

void PolynomialSurfaceTest::constructor()
{
  SurfData sd(fullPath("oneDimQuadratic.txt"));
  PolynomialSurface ps(&sd, 2);
  CPPUNIT_ASSERT(ps.scaler == 0);
  CPPUNIT_ASSERT(ps.xsize == 1);
  CPPUNIT_ASSERT(!ps.builtOK);
  CPPUNIT_ASSERT(!ps.dataModified);
  CPPUNIT_ASSERT(!ps.hasOriginalData());
  CPPUNIT_ASSERT(ps.excludedPoints.size() == 0);
  CPPUNIT_ASSERT(sd == *(ps.sd));
  CPPUNIT_ASSERT(ps.responseIndex == 0);
  CPPUNIT_ASSERT(ps.name == string("Polynomial"));
  CPPUNIT_ASSERT(ps.order == 2);
  CPPUNIT_ASSERT(ps.coefficients.size() == 0);
  CPPUNIT_ASSERT(ps.digits.size() == 2);
  CPPUNIT_ASSERT(ps.termIndex == 0);
}

void PolynomialSurfaceTest::constructorNullData()
{
  PolynomialSurface ps(NULL, 3);
  CPPUNIT_ASSERT(ps.scaler == 0);
  CPPUNIT_ASSERT(ps.xsize == 0);
  CPPUNIT_ASSERT(!ps.builtOK);
  CPPUNIT_ASSERT(!ps.dataModified);
  CPPUNIT_ASSERT(!ps.hasOriginalData());
  CPPUNIT_ASSERT(ps.excludedPoints.size() == 0);
  CPPUNIT_ASSERT(ps.sd == 0);
  CPPUNIT_ASSERT(ps.responseIndex == 0);
  CPPUNIT_ASSERT(ps.name == string("Polynomial"));
  CPPUNIT_ASSERT(ps.order == 3);
  CPPUNIT_ASSERT(ps.coefficients.size() == 0);
  CPPUNIT_ASSERT(ps.digits.size() == 3);
  CPPUNIT_ASSERT(ps.termIndex == 0);
}

void PolynomialSurfaceTest::constructorZeroOrder()
{
  SurfData sd(fullPath("oneDimQuadratic.txt"));
  PolynomialSurface ps(&sd, 0);
  CPPUNIT_ASSERT(ps.scaler == 0);
  CPPUNIT_ASSERT(ps.xsize == 1);
  CPPUNIT_ASSERT(!ps.builtOK);
  CPPUNIT_ASSERT(!ps.dataModified);
  CPPUNIT_ASSERT(!ps.hasOriginalData());
  CPPUNIT_ASSERT(ps.excludedPoints.size() == 0);
  CPPUNIT_ASSERT(sd == *(ps.sd));
  CPPUNIT_ASSERT(ps.responseIndex == 0);
  CPPUNIT_ASSERT(ps.name == string("Polynomial"));
  CPPUNIT_ASSERT(ps.order == 0);
  CPPUNIT_ASSERT(ps.coefficients.size() == 0);
  CPPUNIT_ASSERT(ps.digits.size() == 0);
  CPPUNIT_ASSERT(ps.termIndex == 0);
}

void PolynomialSurfaceTest::constructorFromTextFile()
{
  SurfData sd(fullPath("oneDimQuadratic.txt"));
  PolynomialSurface ps(fullPath("oneDQpoly2.txt"));
  CPPUNIT_ASSERT(ps.scaler == 0);
  CPPUNIT_ASSERT(ps.xsize == 1);
  CPPUNIT_ASSERT(ps.builtOK);
  CPPUNIT_ASSERT(!ps.dataModified);
  CPPUNIT_ASSERT(ps.hasOriginalData());
  CPPUNIT_ASSERT(ps.excludedPoints.size() == 0);
  CPPUNIT_ASSERT(sd == *(ps.sd));
  CPPUNIT_ASSERT(ps.responseIndex == 0);
  CPPUNIT_ASSERT(ps.name == string("Polynomial"));
  CPPUNIT_ASSERT(ps.order == 2);
  CPPUNIT_ASSERT(ps.coefficients.size() == 3);
  CPPUNIT_ASSERT(ps.digits.size() == 2);
}

void PolynomialSurfaceTest::constructorFromBinaryFile()
{
  SurfData sd(fullPath("oneDimQuadratic.txt"));
  PolynomialSurface ps(fullPath("oneDQpoly2.srf"));
  CPPUNIT_ASSERT(ps.scaler == 0);
  CPPUNIT_ASSERT(ps.xsize == 1);
  CPPUNIT_ASSERT(ps.builtOK);
  CPPUNIT_ASSERT(!ps.dataModified);
  CPPUNIT_ASSERT(ps.hasOriginalData());
  CPPUNIT_ASSERT(ps.excludedPoints.size() == 0);
  CPPUNIT_ASSERT(sd == *(ps.sd));
  CPPUNIT_ASSERT(ps.responseIndex == 0);
  CPPUNIT_ASSERT(ps.name == string("Polynomial"));
  CPPUNIT_ASSERT(ps.order == 2);
  CPPUNIT_ASSERT(ps.coefficients.size() == 3);
  CPPUNIT_ASSERT(ps.digits.size() == 2);
}

void PolynomialSurfaceTest::constructorExplicit()
{
  vector<double> coefficients(3);
  coefficients[0] = 0.0;
  coefficients[1] = 0.0;
  coefficients[2] = 1.0;
  PolynomialSurface ps(1,2,coefficients);
  CPPUNIT_ASSERT(ps.scaler == 0);
  CPPUNIT_ASSERT(ps.xsize == 1);
  CPPUNIT_ASSERT(ps.builtOK);
  CPPUNIT_ASSERT(!ps.dataModified);
  CPPUNIT_ASSERT(!ps.hasOriginalData());
  CPPUNIT_ASSERT(ps.excludedPoints.size() == 0);
  CPPUNIT_ASSERT(ps.sd == 0);
  CPPUNIT_ASSERT(ps.responseIndex == 0);
  CPPUNIT_ASSERT(ps.name == string("Polynomial"));
  CPPUNIT_ASSERT(ps.order == 2);
  CPPUNIT_ASSERT(ps.coefficients.size() == 3);
  // While it is not normal to use the == operator on doubles, these values should
  // be exactly what they were explicitly set to
  CPPUNIT_ASSERT(ps.coefficients[0] == 0.0);
  CPPUNIT_ASSERT(ps.coefficients[1] == 0.0);
  CPPUNIT_ASSERT(ps.coefficients[2] == 1.0);
  CPPUNIT_ASSERT(ps.digits.size() == 2);
}

void PolynomialSurfaceTest::constructorCopy()
{
  vector<double> coefficients(3);
  coefficients[0] = 0.0;
  coefficients[1] = 0.0;
  coefficients[2] = 1.0;
  PolynomialSurface ps2(1,2,coefficients);
  PolynomialSurface ps(ps2);
  CPPUNIT_ASSERT(ps.scaler == 0);
  CPPUNIT_ASSERT(ps.xsize == 1);
  CPPUNIT_ASSERT(ps.builtOK);
  CPPUNIT_ASSERT(!ps.dataModified);
  CPPUNIT_ASSERT(!ps.hasOriginalData());
  CPPUNIT_ASSERT(ps.excludedPoints.size() == 0);
  CPPUNIT_ASSERT(ps.sd == 0);
  CPPUNIT_ASSERT(ps.responseIndex == 0);
  CPPUNIT_ASSERT(ps.name == string("Polynomial"));
  CPPUNIT_ASSERT(ps.order == 2);
  CPPUNIT_ASSERT(ps.coefficients.size() == 3);
  // While it is not normal to use the == operator on doubles, these values should
  // be exactly what they were explicitly set to
  CPPUNIT_ASSERT(ps.coefficients[0] == 0.0);
  CPPUNIT_ASSERT(ps.coefficients[1] == 0.0);
  CPPUNIT_ASSERT(ps.coefficients[2] == 1.0);
  CPPUNIT_ASSERT(ps.digits.size() == 2);
}

void PolynomialSurfaceTest::makeSimilarWithNewData()
{
  SurfData sd(fullPath("oneDimQuadratic.txt"));
  vector<double> coefficients(3);
  coefficients[0] = 0.0;
  coefficients[1] = 0.0;
  coefficients[2] = 1.0;
  PolynomialSurface ps2(1,2,coefficients);
  PolynomialSurface* ps = ps2.makeSimilarWithNewData(&sd);
  CPPUNIT_ASSERT(ps->xsize == 1);
  CPPUNIT_ASSERT(ps->builtOK);
  CPPUNIT_ASSERT(!ps->dataModified);
  CPPUNIT_ASSERT(ps->hasOriginalData());
  CPPUNIT_ASSERT(ps->excludedPoints.size() == 0);
  CPPUNIT_ASSERT(*(ps->sd) == sd);
  CPPUNIT_ASSERT(ps->responseIndex == 0);
  CPPUNIT_ASSERT(ps->name == string("Polynomial"));
  CPPUNIT_ASSERT(ps->order == 2);
  CPPUNIT_ASSERT(ps->coefficients.size() == 3);
  // While it is not normal to use the == operator on doubles, these values should
  // be exactly what they were explicitly set to
  CPPUNIT_ASSERT(ps->coefficients[0] == 0.0);
  CPPUNIT_ASSERT(ps->coefficients[1] == 0.0);
  CPPUNIT_ASSERT(ps->coefficients[2] == 1.0);
  CPPUNIT_ASSERT(ps->digits.size() == 2);
  delete ps;
  ps = 0;
}

void PolynomialSurfaceTest::surfaceName()
{
  vector<double> coefficients(3);
  coefficients[0] = 0.0;
  coefficients[1] = 0.0;
  coefficients[2] = 1.0;
  PolynomialSurface ps2(1,2,coefficients);
  CPPUNIT_ASSERT(ps2.surfaceName() == string("Polynomial"));
}

//void PolynomialSurfaceTest::minPointsRequiredStatic()
//{
//  CPPUNIT_ASSERT(PolynomialSurface::minPointsRequired(1,1) == 2);
//  CPPUNIT_ASSERT(PolynomialSurface::minPointsRequired(1,1) == 2);
//  CPPUNIT_ASSERT(PolynomialSurface::minPointsRequired(1,1) == 2);
//  CPPUNIT_ASSERT(PolynomialSurface::minPointsRequired(1,1) == 2);
//  CPPUNIT_ASSERT(PolynomialSurface::minPointsRequired(1,1) == 2);
//  CPPUNIT_ASSERT(PolynomialSurface::minPointsRequired(1,1) == 2);
//}

void PolynomialSurfaceTest::resetTermCounter()
{
  vector<double> coefficients(3);
  coefficients[0] = 0.0;
  coefficients[1] = 0.0;
  coefficients[2] = 1.0;
  PolynomialSurface ps2(1,2,coefficients);
  ps2.resetTermCounter();
  ps2.nextTerm();
  ps2.nextTerm();
  ps2.nextTerm();
  ps2.nextTerm();
  ps2.resetTermCounter();
  CPPUNIT_ASSERT(ps2.digits[0] == 0);
  CPPUNIT_ASSERT(ps2.digits[1] == 0);
  CPPUNIT_ASSERT(ps2.digits.size() == 2);
}

//void PolynomialSurfaceTest::fact()
//{
//  CPPUNIT_ASSERT( PolynomialSurface::fact(0) == 1 );
//  CPPUNIT_ASSERT( PolynomialSurface::fact(1) == 1 );
//  CPPUNIT_ASSERT( PolynomialSurface::fact(2) == 2 );
//  CPPUNIT_ASSERT( PolynomialSurface::fact(3) == 6 );
//  CPPUNIT_ASSERT( PolynomialSurface::fact(4) == 24 );
//  CPPUNIT_ASSERT( PolynomialSurface::fact(12) == 479001600 );
//}
//
//void PolynomialSurfaceTest::nChooseR()
//{
//  CPPUNIT_ASSERT( PolynomialSurface::nChooseR(10,6) == 210); 
//}

void PolynomialSurfaceTest::minPointsRequiredStatic()
{
  CPPUNIT_ASSERT( PolynomialSurface::minPointsRequired(10,5) == 3003);
  CPPUNIT_ASSERT( PolynomialSurface::minPointsRequired(3,4) == 35);
  CPPUNIT_ASSERT( PolynomialSurface::minPointsRequired(2,2) == 6);
  CPPUNIT_ASSERT( PolynomialSurface::minPointsRequired(15,10) == 3268760);
  CPPUNIT_ASSERT( PolynomialSurface::minPointsRequired(1,0) == 1);
}

void PolynomialSurfaceTest::minPointsRequiredNonStatic()
{
  vector<double> coefficients;
  PolynomialSurface ps2(4,5,coefficients);
  CPPUNIT_ASSERT(ps2.minPointsRequired() == 126);
}

void PolynomialSurfaceTest::minPointsRequiredNullException()
{
  PolynomialSurface::minPointsRequired(0, 2);
}

void PolynomialSurfaceTest::evaluate()
{
  vector<double> coefficients(6);
  coefficients[0] = 10.0;
  coefficients[1] = 4.0;
  coefficients[2] = 0.0;
  coefficients[3] = 0.0;
  coefficients[4] = -1.0;
  coefficients[5] = 2.0;
  PolynomialSurface ps2(2,2,coefficients);
  vector<double> x(2);
  x[0] = 2.0;
  x[1] = -1.0;
  CPPUNIT_ASSERT(matches(ps2.evaluate(x),22.0)); 
}

void PolynomialSurfaceTest::build()
{
  SurfData sd(fullPath("oneDimQuadratic.txt"));
  PolynomialSurface ps(&sd, 2);
  ps.createModel();
  ps.write(fullPath("oneDQpoly2.txt"));
  CPPUNIT_ASSERT(ps.coefficients.size() == 3);
  CPPUNIT_ASSERT(matches(ps.coefficients[0],0.0));
  CPPUNIT_ASSERT(matches(ps.coefficients[1],0.0));
  CPPUNIT_ASSERT(matches(ps.coefficients[2],1.0));
}

void PolynomialSurfaceTest::computeTerm()
{
  SurfData sd(fullPath("oneDimQuadratic.txt"));
  PolynomialSurface ps(&sd, 2);
  ps.createModel();
  ps.resetTermCounter();
  vector <double> x(1);
  x[0] = 4.0;
  CPPUNIT_ASSERT( matches(ps.computeTerm(x),1.0));  
  ps.nextTerm();
  CPPUNIT_ASSERT( matches(ps.computeTerm(x),4.0));  
  ps.nextTerm();
  CPPUNIT_ASSERT( matches(ps.computeTerm(x),16.0));  

}

void PolynomialSurfaceTest::computeTermException()
{
  SurfData sd(fullPath("oneDimQuadratic.txt"));
  PolynomialSurface ps(&sd, 2);
  ps.createModel();
  ps.resetTermCounter();
  vector <double> x(1);
  x[0] = 4.0;
  ps.nextTerm();
  ps.digits[0] = 2;
  ps.computeTerm(x);  

}
void PolynomialSurfaceTest::nextTerm()
{
  vector<double> coefficients;
  PolynomialSurface ps(4,5,coefficients);
  ps.resetTermCounter();
  // check a term in the middle
  unsigned i;
  for (i = 0; i < 9; i++) {
    ps.nextTerm();
  }
  CPPUNIT_ASSERT( ps.digits.size() == 5 );
  CPPUNIT_ASSERT( ps.digits[0] == 2 );
  CPPUNIT_ASSERT( ps.digits[1] == 2 );
  CPPUNIT_ASSERT( ps.digits[2] == 0 );
  CPPUNIT_ASSERT( ps.digits[3] == 0 );
  CPPUNIT_ASSERT( ps.digits[4] == 0 );
  CPPUNIT_ASSERT( ps.termIndex == 9 );
  // check the last term
  for (i = 0; i < 116; i++) {
    ps.nextTerm();
  }
  CPPUNIT_ASSERT( ps.digits[0] == 4 );
  CPPUNIT_ASSERT( ps.digits[1] == 4 );
  CPPUNIT_ASSERT( ps.digits[2] == 4 );
  CPPUNIT_ASSERT( ps.digits[3] == 4 );
  CPPUNIT_ASSERT( ps.digits[4] == 4 );
  CPPUNIT_ASSERT_EQUAL( ps.termIndex, (unsigned) 125 );
  // Check to make sure termIndex stops getting incremented once 
  // the last term is reached
  for (i = 0; i < 16; i++) {
    ps.nextTerm();
  }
  CPPUNIT_ASSERT( ps.termIndex == 125 );

}

void PolynomialSurfaceTest::nextTermException()
{
  //PolynomialSurface::minPointsRequired(195,5); 
  vector<double> coefficients;
  PolynomialSurface ps(195, 5, coefficients);
  ps.digits[0] = 154;
  ps.digits[1] = 147;
  ps.digits[2] = 138;
  ps.digits[3] = 137;
  ps.digits[4] = 61;
  ps.termIndex = INT_MAX - 1;
  ps.nextTerm();
  ps.nextTerm();
}

void PolynomialSurfaceTest::derivative_1d_0o()
{
  vector<double> coefficients(1);
  coefficients[0] = 3.0;
  PolynomialSurface ps(1,0,coefficients);
  vector<double> gradient(1);
  vector<double> testPoint(1);
  testPoint[0] = 1.0;
  CPPUNIT_ASSERT(matches(ps.getValue(testPoint), 3.0));
  ps.gradient(testPoint,gradient);
  CPPUNIT_ASSERT(matches(gradient[0],0.0));
  testPoint[0] = -1.0;
  ps.gradient(testPoint,gradient);
  CPPUNIT_ASSERT(matches(gradient[0],0.0));
  coefficients[0] = -4.0;
  ps = PolynomialSurface(1,0,coefficients);
  CPPUNIT_ASSERT(matches(ps.getValue(testPoint), -4.0));
  ps.gradient(testPoint,gradient);
  CPPUNIT_ASSERT(matches(gradient[0],0.0));
}

void PolynomialSurfaceTest::derivative_1d_2o()
{
  vector<double> coefficients(3);
  coefficients[0] = 3.0;
  coefficients[1] = 2.0;
  coefficients[2] = 5.0;
  PolynomialSurface ps(1,2,coefficients);
  vector<double> gradient(1);
  vector<double> testPoint(1);
  testPoint[0] = 1.0;
  CPPUNIT_ASSERT(matches(ps.getValue(testPoint), 10.0));
  ps.gradient(testPoint,gradient);
  CPPUNIT_ASSERT(matches(gradient[0],12.0));
  testPoint[0] = -1.0;
  ps.gradient(testPoint,gradient);
  CPPUNIT_ASSERT(matches(gradient[0],-8.0));
  coefficients[0] = -3.0;
  coefficients[1] = -2.0;
  ps = PolynomialSurface(1,2,coefficients);
  ps.gradient(testPoint,gradient);
  CPPUNIT_ASSERT(matches(gradient[0],-12.0));
}

void PolynomialSurfaceTest::derivative_2d_0o()
{
  vector<double> coefficients(1);
  coefficients[0] = 3.0;
  PolynomialSurface ps(2,0,coefficients);
  vector<double> gradient(2);
  vector<double> testPoint(2);
  testPoint[0] = 1.0;
  testPoint[1] = 1.0;
  CPPUNIT_ASSERT(matches(ps.getValue(testPoint), 3.0));
  ps.gradient(testPoint,gradient);
  CPPUNIT_ASSERT(matches(gradient[0],0.0));
  CPPUNIT_ASSERT(matches(gradient[1],0.0));
  testPoint[0] = -1.0;
  testPoint[1] = -1.0;
  ps.gradient(testPoint,gradient);
  CPPUNIT_ASSERT(matches(gradient[0],0.0));
  CPPUNIT_ASSERT(matches(gradient[1],0.0));
  coefficients[0] = -4.0;
  ps = PolynomialSurface(2,0,coefficients);
  CPPUNIT_ASSERT(matches(ps.getValue(testPoint), -4.0));
  ps.gradient(testPoint,gradient);
  CPPUNIT_ASSERT(matches(gradient[0],0.0));
  CPPUNIT_ASSERT(matches(gradient[1],0.0));
}

void PolynomialSurfaceTest::derivative_2d_2o()
{
  // f(x1,x2) = 3 - 2*x1 + x2 + 6*x1*x2 - 4*x2^2
  // df(x1,x2)/dx1 = -2 + 6*x2
  // df(x1,x2)/dx2 = 1+6*x1 - 8*x2
  vector<double> coefficients(6);
  coefficients[0] = 3.0;
  coefficients[1] = -2.0;
  coefficients[2] = 1.0;
  coefficients[3] = 0.0;
  coefficients[4] = 6.0;
  coefficients[5] = -4.0;
  PolynomialSurface ps(2,2,coefficients);
  vector<double> gradient(2);
  vector<double> testPoint(2);
  testPoint[0] = 1.0;
  testPoint[1] = 1.0;
  CPPUNIT_ASSERT(matches(ps.getValue(testPoint), 4.0));
  ps.gradient(testPoint,gradient);
  CPPUNIT_ASSERT(matches(gradient[0],4.0));
  CPPUNIT_ASSERT(matches(gradient[1],-1.0));
  testPoint[0] = -1.0;
  testPoint[1] = -1.0;
  CPPUNIT_ASSERT(matches(ps.getValue(testPoint), 6.0));
  ps.gradient(testPoint,gradient);
  CPPUNIT_ASSERT(matches(gradient[0],-8.0));
  CPPUNIT_ASSERT(matches(gradient[1],3.0));
  coefficients[0] = -3.0;
  coefficients[1] = 2.0;
  ps = PolynomialSurface(2,2,coefficients);
  ps.gradient(testPoint,gradient);
  CPPUNIT_ASSERT(matches(gradient[0],-4.0));
  CPPUNIT_ASSERT(matches(gradient[1],3.0));
}

void PolynomialSurfaceTest::hessian_2d_2o()
{
  vector<double> coefficients(6);
  coefficients[0] = 3.0;
  coefficients[1] = -2.0;
  coefficients[2] = 1.0;
  coefficients[3] = 0.0;
  coefficients[4] = 6.0;
  coefficients[5] = -4.0;
  PolynomialSurface ps(2,2,coefficients);
  vector<double> x(2);
  SurfpackMatrix<double> sm(2,2);
  x[0] = 0.0;
  x[1] = 0.0;
  ps.hessian(x,sm);
  CPPUNIT_ASSERT(matches(sm[0][0],  0.0));
  CPPUNIT_ASSERT(matches(sm[0][1],  6.0));
  CPPUNIT_ASSERT(matches(sm[1][0],  6.0));
  CPPUNIT_ASSERT(matches(sm[1][1], -8.0));
  // Change point to anything else; hessian should be constant
  x[0] = 7.2;
  x[1] = -1.6;
  ps.hessian(x,sm);
  CPPUNIT_ASSERT(matches(sm[0][0],  0.0));
  CPPUNIT_ASSERT(matches(sm[0][1],  6.0));
  CPPUNIT_ASSERT(matches(sm[1][0],  6.0));
  CPPUNIT_ASSERT(matches(sm[1][1], -8.0));
}
void PolynomialSurfaceTest::leastSquares()
{
  // f(x1,x2) = 3 - 2*x1 + x2 + 6*x1*x2 - 4*x2^2
  // df(x1,x2)/dx1 = -2 + 6*x2
  // df(x1,x2)/dx2 = 1+6*x1 - 8*x2
  vector<double> coefficients(6);
  coefficients[0] = 3.0;
  coefficients[1] = -2.0;
  coefficients[2] = 1.0;
  coefficients[3] = 0.0;
  coefficients[4] = 6.0;
  coefficients[5] = -4.0;
  PolynomialSurface ps(2,2,coefficients);
  vector<double> x(2);
  SurfData sd;
  double f;
  x[0] =  0.0; x[1] =  0.0; f =   3.0; sd.addPoint(SurfPoint(x,f));
  x[0] =  1.0; x[1] =  1.0; f =   4.0; sd.addPoint(SurfPoint(x,f));
  x[0] = -1.0; x[1] =  1.0; f =  -4.0; sd.addPoint(SurfPoint(x,f));
  x[0] =  1.0; x[1] = -1.0; f = -10.0; sd.addPoint(SurfPoint(x,f));
  x[0] = -1.0; x[1] = -1.0; f =   6.0; sd.addPoint(SurfPoint(x,f));
  x[0] =  2.0; x[1] =  3.0; f =   2.0; sd.addPoint(SurfPoint(x,f));
  CPPUNIT_ASSERT( matches(ps.getValue(sd[0]),sd.getResponse(0)));
  CPPUNIT_ASSERT( matches(ps.getValue(sd[1]),sd.getResponse(1)));
  CPPUNIT_ASSERT( matches(ps.getValue(sd[2]),sd.getResponse(2)));
  CPPUNIT_ASSERT( matches(ps.getValue(sd[3]),sd.getResponse(3)));
  CPPUNIT_ASSERT( matches(ps.getValue(sd[4]),sd.getResponse(4)));
  CPPUNIT_ASSERT( matches(ps.getValue(sd[5]),sd.getResponse(5)));

  // Now create a surface from the data and see if everything matches
  PolynomialSurface ps2(&sd,2);
  ps2.createModel();
  CPPUNIT_ASSERT( ps2.coefficients.size() == 6 );
  CPPUNIT_ASSERT( matches(ps2.coefficients[0],coefficients[0]));
  CPPUNIT_ASSERT( matches(ps2.coefficients[1],coefficients[1]));
  CPPUNIT_ASSERT( matches(ps2.coefficients[2],coefficients[2]));
  CPPUNIT_ASSERT( matches(ps2.coefficients[3],coefficients[3]));
  CPPUNIT_ASSERT( matches(ps2.coefficients[4],coefficients[4]));
  CPPUNIT_ASSERT( matches(ps2.coefficients[5],coefficients[5]));
}

void PolynomialSurfaceTest::leastSquaresWithConstraints()
{
  // f(x1,x2) = 3 - 2*x1 + x2 + 6*x1*x2 - 4*x2^2
  // df(x1,x2)/dx1 = -2 + 6*x2
  // df(x1,x2)/dx2 = 1+6*x1 - 8*x2
  vector<double> coefficients(6);
  coefficients[0] = 3.0;
  coefficients[1] = -2.0;
  coefficients[2] = 1.0;
  coefficients[3] = 0.0;
  coefficients[4] = 6.0;
  coefficients[5] = -4.0;
  PolynomialSurface ps(2,2,coefficients);
  vector<double> x(2);
  vector<double> gradient(2);
  SurfData sd;
  double f;
  x[0] =  4.0; x[1] =  5.0; f =  20.0; sd.addPoint(SurfPoint(x,f));
  x[0] =  0.0; x[1] =  0.0; f =   3.0; sd.addPoint(SurfPoint(x,f));
  x[0] =  1.0; x[1] =  1.0; f =   4.0; sd.addPoint(SurfPoint(x,f));
  x[0] = -1.0; x[1] =  1.0; f =  -4.0; sd.addPoint(SurfPoint(x,f));
  //x[0] = -1.0; x[1] = -1.0; f =   6.0; sd.addPoint(SurfPoint(x,f));
  //x[0] =  2.0; x[1] =  3.0; f =   2.0; sd.addPoint(SurfPoint(x,f));
  //x[0] =  1.0; x[1] = -1.0; f = -10.0; sd.addPoint(SurfPoint(x,f));
  x[0] =  1.0; x[1] = -1.0; f = -10.0; SurfPoint sp(x,f);
  gradient[0] = -8.0;
  gradient[1] = 15.0;
  
  

  // Now create a surface from the data and see if everything matches
  PolynomialSurface ps2(&sd,2);
  ps2.setEqualityConstraints(3,sp,sp.F(),&gradient,0);
  ps2.createModel();
  CPPUNIT_ASSERT( ps2.coefficients.size() == 6 );
  CPPUNIT_ASSERT( matches(coefficients[0],ps2.coefficients[0]));
  CPPUNIT_ASSERT( matches(coefficients[1],ps2.coefficients[1]));
  CPPUNIT_ASSERT( matches(coefficients[2],ps2.coefficients[2]));
  CPPUNIT_ASSERT( matches(coefficients[3],ps2.coefficients[3]));
  CPPUNIT_ASSERT( matches(coefficients[4],ps2.coefficients[4]));
  CPPUNIT_ASSERT( matches(coefficients[5],ps2.coefficients[5]));

  CPPUNIT_ASSERT( matches(ps.getValue(sd[0]),sd.getResponse(0)));
  //CPPUNIT_ASSERT( matches(ps.getValue(sd[1]),sd.getResponse(1)));
  //CPPUNIT_ASSERT( matches(ps.getValue(sd[2]),sd.getResponse(2)));
}

void PolynomialSurfaceTest::leastSquaresWithHessianConstraints()
{
  // f(x1,x2) = 3 - 2*x1 + x2 + 6*x1*x2 - 4*x2^2
  // df(x1,x2)/dx1 = -2 + 6*x2
  // df(x1,x2)/dx2 = 1+6*x1 - 8*x2
  vector<double> coefficients(6);
  coefficients[0] = 3.0;
  coefficients[1] = -2.0;
  coefficients[2] = 1.0;
  coefficients[3] = 0.0;
  coefficients[4] = 6.0;
  coefficients[5] = -4.0;
  PolynomialSurface ps(2,2,coefficients);
  vector<double> x(2);
  vector<double> gradient(2);
  SurfpackMatrix<double> hessian(2,2);
  hessian[0][0] = 0.0;
  hessian[0][1] = hessian[1][0] = 6.0;
  hessian[1][1] = -8.0;
  SurfData sd;
  double f;
  x[0] =  4.0; x[1] =  5.0; f =  20.0; sd.addPoint(SurfPoint(x,f));
  x[0] =  0.0; x[1] =  0.0; f =   3.0; sd.addPoint(SurfPoint(x,f));
  x[0] =  1.0; x[1] =  1.0; f =   4.0; sd.addPoint(SurfPoint(x,f));
  x[0] = -1.0; x[1] =  1.0; f =  -4.0; sd.addPoint(SurfPoint(x,f));
  //x[0] = -1.0; x[1] = -1.0; f =   6.0; sd.addPoint(SurfPoint(x,f));
  //x[0] =  2.0; x[1] =  3.0; f =   2.0; sd.addPoint(SurfPoint(x,f));
  //x[0] =  1.0; x[1] = -1.0; f = -10.0; sd.addPoint(SurfPoint(x,f));
  x[0] =  1.0; x[1] = -1.0; f = -10.0; SurfPoint sp(x,f);
  gradient[0] = -8.0;
  gradient[1] = 15.0;
  
  

  // Now create a surface from the data and see if everything matches
  PolynomialSurface ps2(&sd,2);
  ps2.setEqualityConstraints(7,sp,sp.F(),&gradient,&hessian);
  // Spot check the values in the matrix 
  CPPUNIT_ASSERT( matches(ps2.eqConLHS[3][3],2.0));
  CPPUNIT_ASSERT( matches(ps2.eqConLHS[4][4],1.0));
  CPPUNIT_ASSERT( matches(ps2.eqConLHS[5][5],2.0));
  CPPUNIT_ASSERT( matches(ps2.eqConLHS[4][3],0.0));
  CPPUNIT_ASSERT( matches(ps2.eqConLHS[5][4],0.0));
  CPPUNIT_ASSERT( matches(ps2.eqConLHS[3][5],0.0));
  CPPUNIT_ASSERT( matches(ps2.eqConRHS[3], 0.0));
  CPPUNIT_ASSERT( matches(ps2.eqConRHS[4], 6.0));
  CPPUNIT_ASSERT( matches(ps2.eqConRHS[5],-8.0));
  ps2.createModel();
  CPPUNIT_ASSERT( ps2.coefficients.size() == 6 );
  CPPUNIT_ASSERT( matches(coefficients[0],ps2.coefficients[0]));
  CPPUNIT_ASSERT( matches(coefficients[1],ps2.coefficients[1]));
  CPPUNIT_ASSERT( matches(coefficients[2],ps2.coefficients[2]));
  CPPUNIT_ASSERT( matches(coefficients[3],ps2.coefficients[3]));
  CPPUNIT_ASSERT( matches(coefficients[4],ps2.coefficients[4]));
  CPPUNIT_ASSERT( matches(coefficients[5],ps2.coefficients[5]));

  CPPUNIT_ASSERT( matches(ps.getValue(sd[0]),sd.getResponse(0)));

  //CPPUNIT_ASSERT( matches(ps.getValue(sd[1]),sd.getResponse(1)));
  //CPPUNIT_ASSERT( matches(ps.getValue(sd[2]),sd.getResponse(2)));
}

void PolynomialSurfaceTest::leastSquares_2d_3o()
{
  // f(x1,x2) = 0 - x1 + 2*x2 - 3*x1^2 + 4*x1*x2 - 5*x2^2 + 6*x1^3 
		//- 7*x1^2*x2 + 8*x1*x2^2 - 9*x2^3 
  // df(x1,x2)/dx1 = -1 -6*x1 + 4*x2 + 18*x1^2 - 14*x1*x2 + 8*x2^2 
  // df(x1,x2)/dx2 = 2 + 4*x1 - 10*x2 - 7*x1^2 + 16*x1*x2 - 27*x2^2 
  // df^2(x1,x2)/dx1^2 = -6 + 36*x1 - 14*x2
  // df^2(x1,x2)/dx2^2 = -10 + 16*x1 - 54*x2
  // df^2(x1,x2)/dx1x2 = 4 - 14*x1 + 16*x2
 
  vector<double> coefficients(10);
  coefficients[0] =  0.0;
  coefficients[1] = -1.0;
  coefficients[2] =  2.0;
  coefficients[3] = -3.0;
  coefficients[4] =  4.0;
  coefficients[5] = -5.0;
  coefficients[6] =  6.0;
  coefficients[7] = -7.0;
  coefficients[8] =  8.0;
  coefficients[9] = -9.0;
  PolynomialSurface ps(2,3,coefficients);
  vector<double> x(2);
  vector<double> gradient(2);
  SurfpackMatrix<double> hessian(2,2);
  x[0] = 0.0;
  x[1] = 0.0;

  // check the value, gradient, and hessian at this point
  CPPUNIT_ASSERT( matches(ps.getValue(x), 0.0));
  ps.gradient(x,gradient);
  CPPUNIT_ASSERT( matches(gradient[0], -1.0));
  CPPUNIT_ASSERT( matches(gradient[1],  2.0));
  ps.hessian(x,hessian);
  CPPUNIT_ASSERT( matches(hessian[0][0], -6.0));
  CPPUNIT_ASSERT( matches(hessian[0][1],  4.0));
  CPPUNIT_ASSERT( matches(hessian[1][0],  4.0));
  CPPUNIT_ASSERT( matches(hessian[1][1], -10.0));

  // now pick another point and do all the same checks
  x[0] = 10.0;
  x[1] = -7.5;

  // check the value, gradient, and hessian at this point
  CPPUNIT_ASSERT( matches(ps.getValue(x), 18640.625));
  ps.gradient(x,gradient);
  CPPUNIT_ASSERT( matches(gradient[0], 3209.0));
  CPPUNIT_ASSERT( matches(gradient[1],  -3301.75));
  ps.hessian(x,hessian);
  CPPUNIT_ASSERT( matches(hessian[0][0], 459.0));
  CPPUNIT_ASSERT( matches(hessian[0][1], -256.0));
  CPPUNIT_ASSERT( matches(hessian[1][0], -256.0));
  CPPUNIT_ASSERT( matches(hessian[1][1], 555.0));
  
  // Evaluate this surface at 10 points and use those points to
  // create another surface-- in theory, with the same coefficients
  SurfData sd;
  double f;
  x[0] =  6.0; x[1] = -1.0; f =  ps.getValue(x); sd.addPoint(SurfPoint(x,f));
  x[0] = 14.0; x[1] = 35.0; f =  ps.getValue(x); sd.addPoint(SurfPoint(x,f));
  x[0] =  2.0; x[1] =  8.0; f =  ps.getValue(x); sd.addPoint(SurfPoint(x,f));
  x[0] = 21.0; x[1] = -5.0; f =  ps.getValue(x); sd.addPoint(SurfPoint(x,f));
  x[0] = 24.0; x[1] = 15.0; f =  ps.getValue(x); sd.addPoint(SurfPoint(x,f));
  x[0] =  4.0; x[1] = -5.0; f =  ps.getValue(x); sd.addPoint(SurfPoint(x,f));
  x[0] =  0.0; x[1] = 25.0; f =  ps.getValue(x); sd.addPoint(SurfPoint(x,f));
  x[0] = -4.0; x[1] =  9.0; f =  ps.getValue(x); sd.addPoint(SurfPoint(x,f));
  x[0] =  8.0; x[1] =  3.0; f =  ps.getValue(x); sd.addPoint(SurfPoint(x,f));
  x[0] = -9.0; x[1] = -2.0; f =  ps.getValue(x); sd.addPoint(SurfPoint(x,f));
  x[0] =  3.0; x[1] =  1.0; f =  ps.getValue(x); sd.addPoint(SurfPoint(x,f));
  PolynomialSurface ps2(&sd,3);
  ps2.createModel();
  CPPUNIT_ASSERT( matches( ps2.coefficients[0],coefficients[0] ));
  CPPUNIT_ASSERT( matches( ps2.coefficients[1],coefficients[1] ));
  CPPUNIT_ASSERT( matches( ps2.coefficients[2],coefficients[2] ));
  CPPUNIT_ASSERT( matches( ps2.coefficients[3],coefficients[3] ));
  CPPUNIT_ASSERT( matches( ps2.coefficients[4],coefficients[4] ));
  CPPUNIT_ASSERT( matches( ps2.coefficients[5],coefficients[5] ));
  CPPUNIT_ASSERT( matches( ps2.coefficients[6],coefficients[6] ));
  CPPUNIT_ASSERT( matches( ps2.coefficients[7],coefficients[7] ));
  CPPUNIT_ASSERT( matches( ps2.coefficients[8],coefficients[8] ));
  CPPUNIT_ASSERT( matches( ps2.coefficients[9],coefficients[9] ));
  
  // Okay, now try it using a few points, plus deriv and hessian info
  SurfData sd2;
  x[0] =  6.0; x[1] = -1.0; f =  ps.getValue(x); sd2.addPoint(SurfPoint(x,f));
  x[0] = 14.0; x[1] = 35.0; f =  ps.getValue(x); sd2.addPoint(SurfPoint(x,f));
  x[0] =  2.0; x[1] =  8.0; f =  ps.getValue(x); sd2.addPoint(SurfPoint(x,f));
  x[0] = 21.0; x[1] = -5.0; f =  ps.getValue(x); sd2.addPoint(SurfPoint(x,f));
  //x[0] = 24.0; x[1] = 15.0; f =  ps.getValue(x); sd2.addPoint(SurfPoint(x,f));
  //x[0] =  4.0; x[1] = -5.0; f =  ps.getValue(x); sd2.addPoint(SurfPoint(x,f));
  //x[0] =  0.0; x[1] = 25.0; f =  ps.getValue(x); sd2.addPoint(SurfPoint(x,f));

  // Now here is the point the model must fit exactly (resp,grad,hess)
  x[0] = 3.0; x[1] = 1.0; f =  ps.getValue(x); 
  CPPUNIT_ASSERT(matches(f,93.0));
  SurfPoint sp2(x,f);
  ps.gradient(x,gradient);
  CPPUNIT_ASSERT(matches(gradient[0],113.0));
  CPPUNIT_ASSERT(matches(gradient[1],-38.0));
  ps.hessian(x,hessian);
  CPPUNIT_ASSERT(matches(hessian[0][0],88.0));
  CPPUNIT_ASSERT(matches(hessian[0][1],-22.0));
  CPPUNIT_ASSERT(matches(hessian[1][0],-22.0));
  CPPUNIT_ASSERT(matches(hessian[1][1],-16.0));

  PolynomialSurface ps3(&sd2,3);
  ps3.setEqualityConstraints(7,sp2,f,&gradient,&hessian);
  // Spot check the values in the matrix 
  CPPUNIT_ASSERT( matches(ps3.eqConLHS[0][3],9.0));
  CPPUNIT_ASSERT( matches(ps3.eqConLHS[0][6],27.0));
  CPPUNIT_ASSERT( matches(ps3.eqConLHS[1][6],27.0));
  CPPUNIT_ASSERT( matches(ps3.eqConLHS[1][7],6.0));
  CPPUNIT_ASSERT( matches(ps3.eqConLHS[2][4],3.0));
  CPPUNIT_ASSERT( matches(ps3.eqConLHS[2][9],3.0));
  CPPUNIT_ASSERT( matches(ps3.eqConLHS[3][5],0.0));
  CPPUNIT_ASSERT( matches(ps3.eqConLHS[3][6],18.0));
  CPPUNIT_ASSERT( matches(ps3.eqConLHS[4][3],0.0));
  CPPUNIT_ASSERT( matches(ps3.eqConLHS[4][7],6.0));
  CPPUNIT_ASSERT( matches(ps3.eqConLHS[5][8],6.0));
  CPPUNIT_ASSERT( matches(ps3.eqConLHS[5][9],6.0));
  // Check right-hand sides of constraints
  CPPUNIT_ASSERT( matches(ps3.eqConRHS[0],93.0));
  CPPUNIT_ASSERT( matches(ps3.eqConRHS[1],113.0));
  CPPUNIT_ASSERT( matches(ps3.eqConRHS[2],-38.0));
  CPPUNIT_ASSERT( matches(ps3.eqConRHS[3],88.0));
  CPPUNIT_ASSERT( matches(ps3.eqConRHS[4],-22.0));
  CPPUNIT_ASSERT( matches(ps3.eqConRHS[5],-16.0));
  ps3.createModel();
  CPPUNIT_ASSERT( matches( ps3.coefficients[0],coefficients[0] ));
  CPPUNIT_ASSERT( matches( ps3.coefficients[1],coefficients[1] ));
  CPPUNIT_ASSERT( matches( ps3.coefficients[2],coefficients[2] ));
  CPPUNIT_ASSERT( matches( ps3.coefficients[3],coefficients[3] ));
  CPPUNIT_ASSERT( matches( ps3.coefficients[4],coefficients[4] ));
  CPPUNIT_ASSERT( matches( ps3.coefficients[5],coefficients[5] ));
  CPPUNIT_ASSERT( matches( ps3.coefficients[6],coefficients[6] ));
  CPPUNIT_ASSERT( matches( ps3.coefficients[7],coefficients[7] ));
  CPPUNIT_ASSERT( matches( ps3.coefficients[8],coefficients[8] ));
  CPPUNIT_ASSERT( matches( ps3.coefficients[9],coefficients[9] ));
  
}

void PolynomialSurfaceTest::io()
{
  SurfData sd(fullPath("oneDimQuadratic.txt"));
  PolynomialSurface ps(fullPath("oneDQpoly2.txt"));
  CPPUNIT_ASSERT(ps.xsize == 1);
  CPPUNIT_ASSERT(ps.builtOK);
  CPPUNIT_ASSERT(!ps.dataModified);
  CPPUNIT_ASSERT(ps.hasOriginalData());
  CPPUNIT_ASSERT(ps.excludedPoints.size() == 0);
  CPPUNIT_ASSERT(sd == *(ps.sd));
  CPPUNIT_ASSERT(ps.responseIndex == 0);
  CPPUNIT_ASSERT(ps.name == string("Polynomial"));
  CPPUNIT_ASSERT(ps.order == 2);
  CPPUNIT_ASSERT(ps.coefficients.size() == 3);
  CPPUNIT_ASSERT(ps.digits.size() == 2);
  ps.write(fullPath("oneDQpoly2copy.txt"));
  PolynomialSurface ps2(fullPath("oneDQpoly2copy.txt"));
  ps2.write(fullPath("oneDQpoly2copy.srf"));
  PolynomialSurface ps3(fullPath("oneDQpoly2copy.srf"));
  CPPUNIT_ASSERT(ps3.xsize == 1);
  CPPUNIT_ASSERT(ps3.builtOK);
  CPPUNIT_ASSERT(!ps3.dataModified);
  CPPUNIT_ASSERT(ps3.hasOriginalData());
  CPPUNIT_ASSERT(ps3.excludedPoints.size() == 0);
  CPPUNIT_ASSERT(sd == *(ps3.sd));
  CPPUNIT_ASSERT(ps3.responseIndex == 0);
  CPPUNIT_ASSERT(ps3.name == string("Polynomial"));
  CPPUNIT_ASSERT(ps3.order == 2);
  CPPUNIT_ASSERT(ps3.coefficients.size() == 3);
  CPPUNIT_ASSERT(ps3.digits.size() == 2);
}

void PolynomialSurfaceTest::printTermLabel()
{
  vector<double> coefficients;
  PolynomialSurface ps(4,5,coefficients);
  ps.resetTermCounter();
  unsigned i;
  for (i = 0; i < 9; i++) {
    ps.nextTerm();
  }
  ostringstream ostr;
  ps.printTermLabel(ostr);
  CPPUNIT_ASSERT_EQUAL(ostr.str(),string("x2^2"));
  ps.nextTerm();
  ostringstream ostr2;
  ps.printTermLabel(ostr2);
  CPPUNIT_ASSERT_EQUAL(ostr2.str(),string("x3*x2"));

}
  
void PolynomialSurfaceTest::printTermComponents()
{
  vector<double> coefficients;
  PolynomialSurface ps(4,5,coefficients);
  ps.resetTermCounter();
  unsigned i;
  for (i = 0; i < 9; i++) {
    ps.nextTerm();
  }
  ostringstream ostr;
  ps.printTermComponents(ostr);
  CPPUNIT_ASSERT_EQUAL(ostr.str(),string(" [2][2][][][] "));
  ps.nextTerm();
  ostringstream ostr2;
  ps.printTermComponents(ostr2);
  CPPUNIT_ASSERT_EQUAL(ostr2.str(),string(" [3][2][][][] "));
  cout << "End PolynomialSurfaceTest" << endl;
}
