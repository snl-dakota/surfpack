#include "SurfPointTest.h"
#include "SurfPoint.h"
#include <vector>

using namespace std;

CPPUNIT_TEST_SUITE_REGISTRATION( SurfPointTest );

void SurfPointTest::setUp()
{
}

void SurfPointTest::tearDown()
{
}

void SurfPointTest::testConstructor()
{
  
  vector<double> pt;
  pt.push_back(1);
  SurfPoint s(pt);
  CPPUNIT_ASSERT_EQUAL(s.x[0], 1.0);
}
