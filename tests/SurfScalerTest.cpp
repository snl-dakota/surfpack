#include "config.h"

#include <vector>
#include <fstream>
#include <iostream>
#include <string>

#include "SurfScalerTest.h"
#include "SurfScaler.h"
#include "SurfPoint.h"
#include "SurfData.h"
#include "unittests.h"

using namespace std;

CPPUNIT_TEST_SUITE_REGISTRATION( SurfScalerTest );

void SurfScalerTest::setUp()
{
}

void SurfScalerTest::tearDown()
{
}

void SurfScalerTest::testConstructor()
{
  SurfScaler* ss = new SurfScaler();
  CPPUNIT_ASSERT_EQUAL((unsigned)0, ss->parameters.size());
  CPPUNIT_ASSERT_EQUAL((unsigned)1, ss->scaledPoint.xSize());
  CPPUNIT_ASSERT_EQUAL((double)0, ss->scaledPoint.X()[0]);
  delete ss;
}

//void SurfScalerTest::testConstructorXSpecified()
//{
//}
//void SurfScalerTest::testConstructorXSpecifiedPlusOneF()
//{
//}
//
//void SurfScalerTest::testConstructorXSpecifiedFVector()
//{
//}
//
//void SurfScalerTest::testConstructorFromIStreamBinary()
//{
//}
//
void SurfScalerTest::testCopyConstructor()
{
  SurfScaler* ss = new SurfScaler();
  SurfScaler* ss2 = new SurfScaler(*ss);
  CPPUNIT_ASSERT_EQUAL((double)0, ss2->scaledPoint.X()[0]);
  CPPUNIT_ASSERT_EQUAL(ss->scaledPoint.xSize(), 
    ss2->scaledPoint.xSize());
  CPPUNIT_ASSERT_EQUAL(ss->parameters.size(),
     ss2->parameters.size());
  CPPUNIT_ASSERT_EQUAL((unsigned)0,
     ss2->parameters.size());
  delete ss;
  delete ss2;

}
//
//void SurfScalerTest::testConstructorBadXSize()
//{
//}

// Overloaded operators
void SurfScalerTest::testOperatorAssignment()
{
}

void SurfScalerTest::testOperatorAssignmentToSelf()
{
}

void SurfScalerTest::testComputeScalingParameters()
{
  double i = 1.0;
  vector<double> pt(1);
  vector<double> resp(1);
  vector<SurfPoint> spts;
  for (unsigned j = 0; j < 10; j++) {
    pt[0] = i;
    resp[0] = i;
    spts.push_back(SurfPoint(pt,resp));
    i += 1.0;
  }
  SurfData sd(spts);
  SurfScaler ss;
  ss.computeScalingParameters(sd);
  CPPUNIT_ASSERT_EQUAL((unsigned)1,ss.parameters.size());
  CPPUNIT_ASSERT(matches(ss.parameters[0].offset,1.0));
  CPPUNIT_ASSERT(matches(ss.parameters[0].divisor,10.0));
}
//void SurfScalerTest::testOperatorEquality()
//{
//  SurfPoint sp(x2, f1);
//  SurfPoint sp2(sp);
//  CPPUNIT_ASSERT(sp == sp2);
//  CPPUNIT_ASSERT(sp.operator==(sp2));
//  CPPUNIT_ASSERT(sp2.operator==(sp));
//}
//
//void SurfScalerTest::testOperatorInequality()
//{
//  SurfPoint sp(x2, f1);
//  SurfPoint sp2(x1, 3.0);
//  CPPUNIT_ASSERT(sp != sp2);
//  CPPUNIT_ASSERT(sp.operator!=(sp2));
//  CPPUNIT_ASSERT(sp2.operator!=(sp));
//}

// Queries
// I/O
