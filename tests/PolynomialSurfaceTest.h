#include "config.h"

#ifndef POLYNOMIAL_SURFACE_TEST_H
#define POLYNOMIAL_SURFACE_TEST_H 

#include <cppunit/extensions/HelperMacros.h>
#include <string>

#include "surfpack.h"
#include "PolynomialSurface.h"

class PolynomialSurfaceTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( PolynomialSurfaceTest );
  CPPUNIT_TEST( constructor );
  CPPUNIT_TEST( constructorNullData );
  CPPUNIT_TEST( constructorZeroOrder);
  CPPUNIT_TEST( constructorFromTextFile);
  CPPUNIT_TEST( constructorFromBinaryFile);
  CPPUNIT_TEST( constructorExplicit );
  CPPUNIT_TEST( constructorCopy );
  CPPUNIT_TEST( makeSimilarWithNewData );
  CPPUNIT_TEST( surfaceName );
  CPPUNIT_TEST( minPointsRequiredStatic );
  CPPUNIT_TEST( resetTermCounter );
//  CPPUNIT_TEST( fact );
//  CPPUNIT_TEST( nChooseR );
  CPPUNIT_TEST( minPointsRequiredStatic );
  CPPUNIT_TEST( minPointsRequiredNonStatic );
  CPPUNIT_TEST_EXCEPTION( minPointsRequiredNullException, std::string);
  CPPUNIT_TEST( evaluate );
  CPPUNIT_TEST( build );
  CPPUNIT_TEST( computeTerm );
  CPPUNIT_TEST_EXCEPTION( computeTermException, std::range_error );
  CPPUNIT_TEST( nextTerm );
  CPPUNIT_TEST_EXCEPTION( nextTermException, std::range_error);
  CPPUNIT_TEST( io );
  CPPUNIT_TEST( printTermLabel );
  CPPUNIT_TEST( printTermComponents );
  CPPUNIT_TEST_SUITE_END();
public:
  void setUp();
  void tearDown();
  void constructor();
  void constructorNullData();
  void constructorZeroOrder();
  void constructorFromTextFile();
  void constructorFromBinaryFile();
  void constructorExplicit();
  void constructorCopy();
  void makeSimilarWithNewData();
  void surfaceName();
  void resetTermCounter(); 
//  void fact();
//  void nChooseR();
  void minPointsRequiredStatic();
  void minPointsRequiredNonStatic();
  void minPointsRequiredNullException();
  void evaluate();
  void build();
  void computeTerm();
  void computeTermException();
  void nextTerm();
  void nextTermException();
  void io();
  void printTermLabel();
  void printTermComponents();
private:
};

#endif
