// These lines were automatically prepended
#ifndef SURFPACK_CONFIG_H
#define SURFPACK_CONFIG_H
#include "config.h"
#endif
// End of prepended lines

#include <vector>
#include <set>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>

#include "SurfPoint.h"
#include "SurfData.h"
#include "Surface.h"
#include "unittests.h"
#include "PolynomialSurface.h"
#include "PolynomialSurfaceTest.h"

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
  CPPUNIT_ASSERT(matches(ps2.evaluate(x),40.0)); 
}

void PolynomialSurfaceTest::build()
{
  SurfData sd(fullPath("oneDimQuadratic.txt"));
  PolynomialSurface ps(&sd, 2);
  ps.createModel();
  ps.write("oneDQpoly2.txt");
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
