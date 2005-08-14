#include "surfpack_config.h"

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
  SurfScaler ss;
  CPPUNIT_ASSERT_EQUAL((unsigned)0,ss.scalers.size());
  CPPUNIT_ASSERT_EQUAL((unsigned)0,ss.responseScalers.size());
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
  CPPUNIT_ASSERT_EQUAL((unsigned)1,ss.scalers.size());
  CPPUNIT_ASSERT(matches(dynamic_cast<NormalizingScaler*>(ss.scalers[0])->offset,1.0));
  CPPUNIT_ASSERT(matches(dynamic_cast<NormalizingScaler*>(ss.scalers[0])->divisor,9.0));
}

void SurfScalerTest::testComputeScalingParametersFour1DPts()
{
  vector<double> x(1);
  vector<double> f(1);
  vector<SurfPoint> spts;
  x[0] = -2; f[0] =  3; spts.push_back(SurfPoint(x,f));
  x[0] =  2; f[0] = -3; spts.push_back(SurfPoint(x,f));
  x[0] = -6; f[0] = -6; spts.push_back(SurfPoint(x,f));
  x[0] =  6; f[0] =  6; spts.push_back(SurfPoint(x,f));
  SurfData sd(spts);
  SurfScaler ss;
  ss.computeScalingParameters(spts);
  CPPUNIT_ASSERT_EQUAL((unsigned)1,ss.scalers.size());
  CPPUNIT_ASSERT(matches(dynamic_cast<NormalizingScaler*>(ss.scalers[0])->offset,-6.0));
  CPPUNIT_ASSERT(matches(dynamic_cast<NormalizingScaler*>(ss.scalers[0])->divisor,12.0));
}

void SurfScalerTest::testComputeScalingParametersFour2DPts()
{
  vector<double> x(2);
  vector<double> f(1);
  vector<SurfPoint> spts;
  x[0] = -6; x[1] =  1; f[0] = -6; spts.push_back(SurfPoint(x,f));
  x[0] =  6; x[1] =  2; f[0] =  6; spts.push_back(SurfPoint(x,f));
  x[0] = -2; x[1] =  3; f[0] =  3; spts.push_back(SurfPoint(x,f));
  x[0] =  2; x[1] =  4; f[0] = -3; spts.push_back(SurfPoint(x,f));
  SurfData sd(spts);
  SurfScaler ss;
  ss.computeScalingParameters(spts);
  CPPUNIT_ASSERT_EQUAL((unsigned)2,ss.scalers.size());
  CPPUNIT_ASSERT(matches(dynamic_cast<NormalizingScaler*>(ss.scalers[0])->offset,-6.0));
  CPPUNIT_ASSERT(matches(dynamic_cast<NormalizingScaler*>(ss.scalers[0])->divisor,12.0));
  CPPUNIT_ASSERT(matches(dynamic_cast<NormalizingScaler*>(ss.scalers[1])->offset,1.0));
  CPPUNIT_ASSERT(matches(dynamic_cast<NormalizingScaler*>(ss.scalers[1])->divisor,3.0));
}

void SurfScalerTest::testScale()
{
  vector<double> x(2);
  vector<double> f(1);
  vector<SurfPoint> spts;
  x[0] = -6; x[1] =  1; f[0] = -6; spts.push_back(SurfPoint(x,f));
  x[0] =  6; x[1] =  2; f[0] =  6; spts.push_back(SurfPoint(x,f));
  x[0] = -2; x[1] =  3; f[0] =  3; spts.push_back(SurfPoint(x,f));
  x[0] =  2; x[1] =  4; f[0] = -3; spts.push_back(SurfPoint(x,f));
  SurfData sd(spts);
  SurfScaler ss;
  sd.setScaler(&ss);
  //ss.computeScalingParameters(spts);
  const SurfPoint& temp = sd[3];
  CPPUNIT_ASSERT(matches(temp.X()[0],.66666));
  CPPUNIT_ASSERT(matches(temp.X()[1],1.0));
  cout << "End SurfScaler Test" << endl;
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
