#include <vector>
#include <set>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>

#include "SurfPoint.h"
#include "SurfData.h"
#include "SurfDataTest.h"

using namespace std;

CPPUNIT_TEST_SUITE_REGISTRATION( SurfDataTest );
const unsigned numPoints = 5;
const unsigned dimPoints = 3;
const unsigned unsignedZero = 0;
const double doubleZero = 0.0;
const int intZero = 0;

void SurfDataTest::setUp()
{
  // surfpoints has numPoints points
  // Each point has 3 dimensions and 3 response values
  for (unsigned i = 0; i < numPoints; i++) {
    vector<double> x;
    vector<double> fx;
    for (unsigned j = i; j < i + 3; j++) {
      x.push_back(static_cast<double>(j));
      fx.push_back(static_cast<double>(j*2));
    }
    surfpoints.push_back(SurfPoint(x,fx));
  }
  sdPtr1 = new SurfData(surfpoints);
  sdPtr2 = 0;
}

void SurfDataTest::tearDown()
{
  delete sdPtr1;
  delete sdPtr2;
}

void SurfDataTest::testConstructor()
{
  SurfData sd(surfpoints);
  CPPUNIT_ASSERT_EQUAL(sd.points, surfpoints);
  CPPUNIT_ASSERT_EQUAL(sd.xsize, dimPoints);
  CPPUNIT_ASSERT_EQUAL(sd.fsize, dimPoints);
  CPPUNIT_ASSERT_EQUAL(sd.mapping.size(), numPoints);
  for (unsigned i = 0; i < numPoints; i++) {
    CPPUNIT_ASSERT_EQUAL(sd.mapping[i], i);
  }
  CPPUNIT_ASSERT(sd.excludedPoints.empty());
  CPPUNIT_ASSERT_EQUAL(sd.defaultIndex, unsignedZero);
  CPPUNIT_ASSERT_EQUAL(sd.xMatrix, static_cast<double*>(0));
  CPPUNIT_ASSERT_EQUAL(sd.yVector, static_cast<double*>(0));
  CPPUNIT_ASSERT(!sd.valid.xMatrix);
  CPPUNIT_ASSERT(!sd.valid.yVector);
}

