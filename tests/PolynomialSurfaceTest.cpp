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

void PolynomialSurfaceTest::constructorTest()
{
  CPPUNIT_ASSERT(1 == 1);
}

