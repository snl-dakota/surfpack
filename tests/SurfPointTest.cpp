// These lines were automatically prepended
#ifndef SURFPACK_CONFIG_H
#define SURFPACK_CONFIG_H
#include "config.h"
#endif
// End of prepended lines

#include <vector>
#include <fstream>
#include <iostream>
#include <string>

#include "SurfPointTest.h"
#include "SurfPoint.h"
#include "unittests.h"

using namespace std;

CPPUNIT_TEST_SUITE_REGISTRATION( SurfPointTest );

void SurfPointTest::setUp()
{
  initialize();
  x1.resize(3);
  x1[0] = 3.0;
  x1[1] = -3.0;
  x1[2] = 0.0;

  x2.resize(1);
  x2[0] = 0.0;

  f1.resize(2);
  f1[0] = 1.0;
  f1[1] = -2.0;

  spPtr = new SurfPoint(x1);
  spPtr2 = new SurfPoint(x2, f1);
}

void SurfPointTest::tearDown()
{
  delete spPtr;
  delete spPtr2;
}

void SurfPointTest::testConstructor()
{
  
  vector<double> pt;
  pt.push_back(1);
  SurfPoint s(pt);
  CPPUNIT_ASSERT_EQUAL(s.x[0], 1.0);
}

void SurfPointTest::testConstructorXSpecified()
{
  SurfPoint sp(x1);
  CPPUNIT_ASSERT_EQUAL(sp.x[0], 3.0);
  CPPUNIT_ASSERT_EQUAL(sp.x[1], -3.0);
  CPPUNIT_ASSERT_EQUAL(sp.x[2], 0.0);
  CPPUNIT_ASSERT_EQUAL(sp.x.size(), static_cast<unsigned>(3));
  

}
void SurfPointTest::testConstructorXSpecifiedPlusOneF()
{
  SurfPoint sp(x2, 2.0);
  CPPUNIT_ASSERT_EQUAL(sp.x[0], 0.0);
  CPPUNIT_ASSERT_EQUAL(sp.f[0], 2.0);
  CPPUNIT_ASSERT_EQUAL(sp.x.size(), static_cast<unsigned>(1));
  CPPUNIT_ASSERT_EQUAL(sp.f.size(), static_cast<unsigned>(1));
}

void SurfPointTest::testConstructorXSpecifiedFVector()
{
  SurfPoint sp(x2, f1);
  CPPUNIT_ASSERT_EQUAL(sp.x[0], 0.0);
  CPPUNIT_ASSERT_EQUAL(sp.f[0], 1.0);
  CPPUNIT_ASSERT_EQUAL(sp.f[1], -2.0);
  CPPUNIT_ASSERT_EQUAL(sp.x.size(), static_cast<unsigned>(1));
  CPPUNIT_ASSERT_EQUAL(sp.f.size(), static_cast<unsigned>(2));

}

void SurfPointTest::testConstructorFromIStreamBinary()
{
  // Create the binary SurfPoint in a file that will be read later
  //vector<double> xs(2);
  //xs[0] = 1.0;
  //xs[1] = 2.0;
  //vector<double> fs(2);
  //fs[0] = 3.0;
  //fs[1] = 4.0;
  //SurfPoint sp2(xs, fs);
  //ofstream outfile(fullPath("point1.sp"), ios::out | ios::binary);
  //sp2.writeBinary(outfile);
  //outfile.close();
 
  const string filename = fullPath("point1.sp"); 
  ifstream infile(filename.c_str());
  SurfPoint sp(2, 2, infile, true);  
  infile.close();
  CPPUNIT_ASSERT_EQUAL(sp.x[0], 1.0);
  CPPUNIT_ASSERT_EQUAL(sp.x[1], 2.0);
  CPPUNIT_ASSERT_EQUAL(sp.f[0], 3.0);
  CPPUNIT_ASSERT_EQUAL(sp.f[1], 4.0);
  CPPUNIT_ASSERT_EQUAL(sp.x.size(), static_cast<unsigned>(2));
  CPPUNIT_ASSERT_EQUAL(sp.f.size(), static_cast<unsigned>(2));
}

void SurfPointTest::testConstructorFromIStreamText()
{
  const string filename = fullPath("point1.txt");
  ifstream infile(filename.c_str());
  SurfPoint sp(2, 2, infile, false);  
  infile.close();
  CPPUNIT_ASSERT_EQUAL(sp.x[0], 1.0);
  CPPUNIT_ASSERT_EQUAL(sp.x[1], 2.0);
  CPPUNIT_ASSERT_EQUAL(sp.f[0], 3.0);
  CPPUNIT_ASSERT_EQUAL(sp.f[1], 4.0);
  CPPUNIT_ASSERT_EQUAL(sp.x.size(), static_cast<unsigned>(2));
  CPPUNIT_ASSERT_EQUAL(sp.f.size(), static_cast<unsigned>(2));
}

void SurfPointTest::testCopyConstructor()
{
  SurfPoint sp2(x2, f1);
  SurfPoint sp(sp2);
  CPPUNIT_ASSERT_EQUAL(sp.x[0], 0.0);
  CPPUNIT_ASSERT_EQUAL(sp.f[0], 1.0);
  CPPUNIT_ASSERT_EQUAL(sp.f[1], -2.0);
  CPPUNIT_ASSERT_EQUAL(sp.x.size(), static_cast<unsigned>(1));
  CPPUNIT_ASSERT_EQUAL(sp.f.size(), static_cast<unsigned>(2));
}

void SurfPointTest::testConstructorBadXSize()
{
  vector<double> x;
  SurfPoint sp(x); // should throw an exception
}

// Overloaded operators
void SurfPointTest::testOperatorAssignment()
{
  SurfPoint sp(x2, f1);
  SurfPoint sp2(x1, 3.0);
  sp2 = sp;
  CPPUNIT_ASSERT_EQUAL(sp.x[0], 0.0);
  CPPUNIT_ASSERT_EQUAL(sp.f[0], 1.0);
  CPPUNIT_ASSERT_EQUAL(sp.f[1], -2.0);
  CPPUNIT_ASSERT_EQUAL(sp.x.size(), static_cast<unsigned>(1));
  CPPUNIT_ASSERT_EQUAL(sp.f.size(), static_cast<unsigned>(2));
}

void SurfPointTest::testOperatorAssignmentToSelf()
{
  SurfPoint sp(x2, f1);
  sp = sp;
  CPPUNIT_ASSERT_EQUAL(sp.x[0], 0.0);
  CPPUNIT_ASSERT_EQUAL(sp.f[0], 1.0);
  CPPUNIT_ASSERT_EQUAL(sp.f[1], -2.0);
  CPPUNIT_ASSERT_EQUAL(sp.x.size(), static_cast<unsigned>(1));
  CPPUNIT_ASSERT_EQUAL(sp.f.size(), static_cast<unsigned>(2));
}

void SurfPointTest::testOperatorEquality()
{
  SurfPoint sp(x2, f1);
  SurfPoint sp2(sp);
  CPPUNIT_ASSERT(sp == sp2);
  CPPUNIT_ASSERT(sp.operator==(sp2));
  CPPUNIT_ASSERT(sp2.operator==(sp));
}

void SurfPointTest::testOperatorInequality()
{
  SurfPoint sp(x2, f1);
  SurfPoint sp2(x1, 3.0);
  CPPUNIT_ASSERT(sp != sp2);
  CPPUNIT_ASSERT(sp.operator!=(sp2));
  CPPUNIT_ASSERT(sp2.operator!=(sp));
}

// Queries
void SurfPointTest::testXSize()
{
  CPPUNIT_ASSERT_EQUAL(spPtr->xSize(), static_cast<unsigned>(3));
  CPPUNIT_ASSERT_EQUAL(spPtr2->xSize(), static_cast<unsigned>(1));
}

void SurfPointTest::testFSize()
{
  CPPUNIT_ASSERT_EQUAL(spPtr->fSize(), static_cast<unsigned>(0));
  CPPUNIT_ASSERT_EQUAL(spPtr2->fSize(), static_cast<unsigned>(2));
}

void SurfPointTest::testX()
{
  vector<double> xvec = spPtr->X();
  CPPUNIT_ASSERT_EQUAL(xvec, x1);
}

void SurfPointTest::testFQuery()
{
  CPPUNIT_ASSERT_EQUAL(spPtr2->F(), 1.0);
  CPPUNIT_ASSERT_EQUAL(spPtr2->F(1), -2.0);
}

void SurfPointTest::testFQueryBadIndex()
{
  // should throw an exception, since spPtr has no response values
  double val = spPtr->F(6); 
}

// Commands
void SurfPointTest::testAddResponse()
{
  SurfPoint sp(*spPtr);
  unsigned newIndex = sp.addResponse();
  CPPUNIT_ASSERT_EQUAL(newIndex, static_cast<unsigned>(0));
  CPPUNIT_ASSERT_EQUAL(sp.fSize(), static_cast<unsigned>(1));
  CPPUNIT_ASSERT_EQUAL(sp.F(), 0.0);
  newIndex = sp.addResponse(4.0);
  CPPUNIT_ASSERT_EQUAL(newIndex, static_cast<unsigned>(1));
  CPPUNIT_ASSERT_EQUAL(sp.fSize(), static_cast<unsigned>(2));
  CPPUNIT_ASSERT_EQUAL(sp.F(1), 4.0);
}

void SurfPointTest::testFAssign()
{
  SurfPoint sp(*spPtr2);
  sp.F(1,3.0);
  CPPUNIT_ASSERT_EQUAL(sp.f[1], 3.0);
}

void SurfPointTest::testFAssignBadIndex()
{
  SurfPoint sp(x1, f1);
  // Should throw an exception, since sp has no response values
  sp.F(3, 1.0); 
}

// I/O
void SurfPointTest::testWriteBinary()
{
  ofstream outfile(fullPath("writePoint.sp").c_str(),ios::out | ios::binary);
  spPtr2->writeBinary(outfile);
  outfile.close();
  ifstream infile(fullPath("writePoint.sp").c_str(),ios::in | ios::binary);
  SurfPoint sp(1, 2, infile, true);
  SurfPoint sp2(x2, f1);
  CPPUNIT_ASSERT(sp == sp2);
}

void SurfPointTest::testWriteText()
{
  ofstream outfile(fullPath("writePoint.txt").c_str(),ios::out);
  spPtr2->writeText(outfile);
  outfile.close();
  ifstream infile(fullPath("writePoint.txt").c_str(),ios::in);
  SurfPoint sp(1, 2, infile, false);
  SurfPoint sp2(x2, f1);
  CPPUNIT_ASSERT(sp == sp2);
}

void SurfPointTest::testReadBinary()
{
  const string filename = fullPath("point2.sp");
  ifstream infile(filename.c_str(), ios::in | ios::binary);
  SurfPoint sp(*spPtr2);
  sp.x[0] = 3.0;
  sp.readBinary(infile);
  CPPUNIT_ASSERT(sp == *spPtr2);
  infile.close();
}

void SurfPointTest::testReadText()
{
  const string filename = fullPath("point2.txt");
  ifstream infile(filename.c_str(), ios::in );
  SurfPoint sp(*spPtr2);
  sp.x[0] = 3.0;
  sp.readText(infile);
  CPPUNIT_ASSERT(sp == *spPtr2);
  infile.close();
}

void SurfPointTest::testStreamInsertion()
{
  // If this test doesn't throw an exception,
  // it is presumed to have worked
  ofstream blackhole("/dev/null",ios::out);
  blackhole << (*spPtr) << endl;
  cout << "End SurfPointTest" << endl;
}

