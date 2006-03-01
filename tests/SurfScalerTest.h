/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */
#include "surfpack_config.h"

#ifndef SURFSCALERTEST_H
#define SURFSCALERTEST_H

#include <cppunit/extensions/HelperMacros.h>
#include "SurfPoint.h"

class SurfScalerTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( SurfScalerTest );
  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST( testCopyConstructor );
  //CPPUNIT_TEST_EXCEPTION( testConstructorBadXSize, SurfPoint::null_point );
  CPPUNIT_TEST( testOperatorAssignment );
  CPPUNIT_TEST( testOperatorAssignmentToSelf );
  CPPUNIT_TEST( testComputeScalingParameters);
  CPPUNIT_TEST( testComputeScalingParametersFour1DPts);
  CPPUNIT_TEST( testComputeScalingParametersFour2DPts);
  CPPUNIT_TEST( testScale );
  //CPPUNIT_TEST( testOperatorEquality );
  //CPPUNIT_TEST( testOperatorInequality );
  CPPUNIT_TEST_SUITE_END();
public:
  void setUp();
  void tearDown();

// Constructors
  void testConstructor();
  void testCopyConstructor();

// Overloaded operators
  void testOperatorAssignment();
  void testOperatorAssignmentToSelf();
  void testComputeScalingParameters();
  void testComputeScalingParametersFour1DPts();
  void testComputeScalingParametersFour2DPts();
  void testScale();
  //void testOperatorEquality();
  //void testOperatorInequality();

// Queries

// Commands

// I/O

private:
  
};


#endif
