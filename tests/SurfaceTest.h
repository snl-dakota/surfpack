// These lines were automatically prepended
#ifndef SURFPACK_CONFIG_H
#define SURFPACK_CONFIG_H
#include "config.h"
#endif
// End of prepended lines

#ifndef SURFACE_TEST_H
#define SURFACE_TEST_H 

#include <cppunit/extensions/HelperMacros.h>
#include "surfpack.h"
#include "PolynomialSurface.h"
#include "Surface.h"

class SurfaceTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( SurfaceTest );
  CPPUNIT_TEST( xSize );
  CPPUNIT_TEST( hasOriginalData );
  CPPUNIT_TEST( acceptableData );
  CPPUNIT_TEST_EXCEPTION( acceptableDataExceptionNull, SurfData::bad_surf_data );
  CPPUNIT_TEST_EXCEPTION( acceptableDataExceptionNotEnough, SurfData::bad_surf_data );
  CPPUNIT_TEST( getValueVector );
  CPPUNIT_TEST( getValueSurfPoint );
  CPPUNIT_TEST( getValueSurfData );
  CPPUNIT_TEST( getValueErrorStructs );
  CPPUNIT_TEST( goodnessOfFit );
  CPPUNIT_TEST_EXCEPTION( goodnessOfFitException, Surface::bad_metric );
  CPPUNIT_TEST( press );
  CPPUNIT_TEST_EXCEPTION( pressExceptionInsufficient, SurfData::bad_surf_data );
  CPPUNIT_TEST( rSquared );
  CPPUNIT_TEST( mse );
  CPPUNIT_TEST( sse );
  CPPUNIT_TEST( mrae );
  CPPUNIT_TEST_EXCEPTION( checkDataException, SurfData::bad_surf_data );
  CPPUNIT_TEST( createModelSurfData );
  CPPUNIT_TEST( recreateModel );
  CPPUNIT_TEST_EXCEPTION( writeNoFile, surfpack::file_open_failure);
  CPPUNIT_TEST( writeNoDataText);
  CPPUNIT_TEST( writeNoDataBinary);
  CPPUNIT_TEST_EXCEPTION( readNoFile, surfpack::file_open_failure );
  CPPUNIT_TEST_EXCEPTION( readBadName, surfpack::io_exception);
  CPPUNIT_TEST_EXCEPTION( badFileExtension, surfpack::io_exception);
  CPPUNIT_TEST( print );
  CPPUNIT_TEST_SUITE_END();
public:
  void setUp();
  void tearDown();
  void xSize();
  void hasOriginalData();
  void acceptableData();
  void acceptableDataExceptionNull();
  void acceptableDataExceptionNotEnough();
  void getValueVector();
  void getValueSurfPoint();
  void getValueSurfData();
  void getValueErrorStructs();
  void goodnessOfFit();
  void goodnessOfFitException();
  void press();
  void pressExceptionInsufficient();
  void rSquared();
  void mse();
  void sse();
  void mrae();
  void checkDataException();
  void createModelSurfData();
  void recreateModel();
  void writeNoFile();
  void writeNoDataText();
  void writeNoDataBinary();
  void readNoFile();
  void readBadName();
  void badFileExtension();
  void print();
private:
  SurfData* surfd;
  PolynomialSurface* polysurf;
};

#endif
