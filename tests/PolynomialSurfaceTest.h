/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#ifdef HAVE_CONFIG_H
#include "surfpack_config.h"
#endif

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
  CPPUNIT_TEST( derivative_1d_0o );
  CPPUNIT_TEST( derivative_1d_2o );
  CPPUNIT_TEST( derivative_2d_0o );
  CPPUNIT_TEST( derivative_2d_2o );
  CPPUNIT_TEST( hessian_2d_2o );
  CPPUNIT_TEST( leastSquares );
  CPPUNIT_TEST( leastSquaresWithConstraints );
  CPPUNIT_TEST( leastSquaresWithHessianConstraints );
  CPPUNIT_TEST( leastSquares_2d_3o);
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
  void derivative_1d_0o();
  void derivative_1d_2o();
  void derivative_2d_0o();
  void derivative_2d_2o();
  void hessian_2d_2o();
  void leastSquares();
  void leastSquaresWithConstraints();
  void leastSquaresWithHessianConstraints();
  void leastSquares_2d_3o();
  void io();
  void printTermLabel();
  void printTermComponents();
private:
};

#endif
