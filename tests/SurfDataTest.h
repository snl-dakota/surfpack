#ifndef SURFDATATEST_H
#define SURFDATATEST_H

#include <cppunit/extensions/HelperMacros.h>
#include "SurfData.h"
#include <vector>

class SurfDataTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( SurfDataTest );
  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST_SUITE_END();
public:
  void setUp();
  void tearDown();
  void testConstructor();
private:
  SurfData* sdPtr1;
  SurfData* sdPtr2;
  std::vector<SurfPoint> surfpoints;
};


#endif