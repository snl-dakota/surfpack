// These lines were automatically prepended
#ifndef SURFPACK_CONFIG_H
#define SURFPACK_CONFIG_H
#include "config.h"
#endif
// End of prepended lines

#ifndef POLYNOMIAL_SURFACE_TEST_H
#define POLYNOMIAL_SURFACE_TEST_H 

#include <cppunit/extensions/HelperMacros.h>
#include "surfpack.h"
#include "PolynomialSurface.h"

class PolynomialSurfaceTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( PolynomialSurfaceTest );
  CPPUNIT_TEST_SUITE_END();
public:
  void setUp();
  void tearDown();
private:
};

#endif
