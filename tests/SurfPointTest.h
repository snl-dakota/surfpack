#ifndef SURFPOINTTEST_H
#define SURFPOINTTEST_H

#include <cppunit/extensions/HelperMacros.h>

class SurfPointTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( SurfPointTest );
  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST_SUITE_END();
public:
  void setUp();
  void tearDown();
  void testConstructor();
};


#endif
