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
#include "unittests.h"
#include "SurfDataTest.h"

using namespace std;

CPPUNIT_TEST_SUITE_REGISTRATION( SurfDataTest );
const unsigned numPoints = 5;
const unsigned dimPoints = 3;
const unsigned numResponses = 3;

void SurfDataTest::setUp()
{
  initialize();
  // surfpoints has numPoints points
  // Each point has 3 dimensions and 3 response values
  surfpoints.reserve(numPoints);
  for (unsigned i = 0; i < numPoints; i++) {
    vector<double> x;
    vector<double> fx;
    for (unsigned j = i; j < i + 3; j++) {
      x.push_back(static_cast<double>(j));
      fx.push_back(static_cast<double>(j*2));
    }
    surfpoints.push_back(SurfPoint(x,fx));
  }
  sdPtr1 = new SurfData(surfpoints);
  sdPtr2 = 0;
  
  for (unsigned i = 0; i < numPoints; i++) {
    skipAllPoints.insert(i);
  }
}

void SurfDataTest::tearDown()
{
  delete sdPtr1;
  delete sdPtr2;
}

void SurfDataTest::testConstructorVectorPoints()
{
  SurfData sd(surfpoints);
  CPPUNIT_ASSERT_EQUAL(sd.points.size(), surfpoints.size());
  for (unsigned i = 0; i < surfpoints.size(); i++) {
    CPPUNIT_ASSERT_EQUAL(*sd.points[i], surfpoints[i]);
  }
  CPPUNIT_ASSERT_EQUAL(sd.xsize, dimPoints);
  CPPUNIT_ASSERT_EQUAL(sd.fsize, dimPoints);
  CPPUNIT_ASSERT_EQUAL(sd.mapping.size(), numPoints);
  for (unsigned i = 0; i < numPoints; i++) {
    CPPUNIT_ASSERT_EQUAL(sd.mapping[i], i);
  }
  CPPUNIT_ASSERT(sd.excludedPoints.empty());
  CPPUNIT_ASSERT_EQUAL(sd.defaultIndex, unsignedZero);
  CPPUNIT_ASSERT_EQUAL(sd.xMatrix, static_cast<double*>(0));
  CPPUNIT_ASSERT_EQUAL(sd.yVector, static_cast<double*>(0));
  CPPUNIT_ASSERT(!sd.valid.xMatrix);
  CPPUNIT_ASSERT(!sd.valid.yVector);
}

void SurfDataTest::testConstructorVectorPointsEmpty()
{
  vector<SurfPoint> noPoints;
  SurfData sd(noPoints);
  CPPUNIT_ASSERT_EQUAL(sd.points.size(), noPoints.size());
  for (unsigned i = 0; i < sd.points.size(); i++) {
    CPPUNIT_ASSERT_EQUAL(*sd.points[i], noPoints[i]);
  }
  CPPUNIT_ASSERT_EQUAL(sd.xsize, unsignedZero);
  CPPUNIT_ASSERT_EQUAL(sd.fsize, unsignedZero);
  CPPUNIT_ASSERT_EQUAL(sd.mapping.size(), unsignedZero);
  CPPUNIT_ASSERT(sd.excludedPoints.empty());
  CPPUNIT_ASSERT_EQUAL(sd.defaultIndex, unsignedZero);
  CPPUNIT_ASSERT_EQUAL(sd.xMatrix, static_cast<double*>(0));
  CPPUNIT_ASSERT_EQUAL(sd.yVector, static_cast<double*>(0));
  CPPUNIT_ASSERT(!sd.valid.xMatrix);
  CPPUNIT_ASSERT(!sd.valid.yVector);
}

void SurfDataTest::testConstructorVectorPointsMismatchXs()
{
  vector<double> x1, x2, f1;
  x1.push_back(1.0);
  x1.push_back(1.5);
  x2 = x1;
  x2.push_back(2.0);
  f1.push_back(3.0);
  vector<SurfPoint> datapoints;
  datapoints.push_back(SurfPoint(x1, f1));
  datapoints.push_back(SurfPoint(x2, f1));
  // should throw exception because x1 and x2 have different dimensionality
  SurfData sd(datapoints);
}

void SurfDataTest::testConstructorVectorPointsMismatchFs()
{
  vector<double> x1, f1, f2;
  x1.push_back(1.0);
  x1.push_back(1.5);
  f1.push_back(3.0);
  vector<SurfPoint> datapoints;
  datapoints.push_back(SurfPoint(x1, f1));
  datapoints.push_back(SurfPoint(x1, f2));
  // should throw exception because x1 and x2 have different dimensionality
  SurfData sd(datapoints);
}

void SurfDataTest::testConstructorFilenameText()
{
  unsigned pointsInFile = 100;
  const string filename = fullPath("rast100.txt");
  SurfData sd(filename);
  CPPUNIT_ASSERT_EQUAL(sd.points.size(), pointsInFile);
  CPPUNIT_ASSERT_EQUAL(sd.xsize, static_cast<unsigned>(2));
  CPPUNIT_ASSERT_EQUAL(sd.fsize, static_cast<unsigned>(1));
  CPPUNIT_ASSERT_EQUAL(sd.mapping.size(), pointsInFile);
  for (unsigned i = 0; i < pointsInFile; i++) {
    CPPUNIT_ASSERT_EQUAL(sd.mapping[i], i);
  }
  CPPUNIT_ASSERT(sd.excludedPoints.empty());
  CPPUNIT_ASSERT_EQUAL(sd.defaultIndex, unsignedZero);
  CPPUNIT_ASSERT_EQUAL(sd.xMatrix, static_cast<double*>(0));
  CPPUNIT_ASSERT_EQUAL(sd.yVector, static_cast<double*>(0));
  CPPUNIT_ASSERT(!sd.valid.xMatrix);
  CPPUNIT_ASSERT(!sd.valid.yVector);
}

void SurfDataTest::testConstructorFilenameBinary()
{
  unsigned pointsInFile = 100;
  const string filename = fullPath("rast100.sd");
  SurfData sd(filename);
  CPPUNIT_ASSERT_EQUAL(sd.points.size(), pointsInFile);
  CPPUNIT_ASSERT_EQUAL(sd.xsize, static_cast<unsigned>(2));
  CPPUNIT_ASSERT_EQUAL(sd.fsize, static_cast<unsigned>(1));
  CPPUNIT_ASSERT_EQUAL(sd.mapping.size(), pointsInFile);
  for (unsigned i = 0; i < pointsInFile; i++) {
    CPPUNIT_ASSERT_EQUAL(sd.mapping[i], i);
  }
  CPPUNIT_ASSERT(sd.excludedPoints.empty());
  CPPUNIT_ASSERT_EQUAL(sd.defaultIndex, unsignedZero);
  CPPUNIT_ASSERT_EQUAL(sd.xMatrix, static_cast<double*>(0));
  CPPUNIT_ASSERT_EQUAL(sd.yVector, static_cast<double*>(0));
  CPPUNIT_ASSERT(!sd.valid.xMatrix);
  CPPUNIT_ASSERT(!sd.valid.yVector);
}

void SurfDataTest::testConstructorIStreamText()
{
  unsigned pointsInFile = 100;
  const string filename = fullPath("rast100.txt");
  ifstream infile(filename.c_str(), ios::in);
  SurfData sd(infile, false);
  infile.close();
  CPPUNIT_ASSERT_EQUAL(sd.points.size(), pointsInFile);
  CPPUNIT_ASSERT_EQUAL(sd.xsize, static_cast<unsigned>(2));
  CPPUNIT_ASSERT_EQUAL(sd.fsize, static_cast<unsigned>(1));
  CPPUNIT_ASSERT_EQUAL(sd.mapping.size(), pointsInFile);
  for (unsigned i = 0; i < pointsInFile; i++) {
    CPPUNIT_ASSERT_EQUAL(sd.mapping[i], i);
  }
  CPPUNIT_ASSERT(sd.excludedPoints.empty());
  CPPUNIT_ASSERT_EQUAL(sd.defaultIndex, unsignedZero);
  CPPUNIT_ASSERT_EQUAL(sd.xMatrix, static_cast<double*>(0));
  CPPUNIT_ASSERT_EQUAL(sd.yVector, static_cast<double*>(0));
  CPPUNIT_ASSERT(!sd.valid.xMatrix);
  CPPUNIT_ASSERT(!sd.valid.yVector);
}

void SurfDataTest::testConstructorIStreamBinary()
{
  unsigned pointsInFile = 100;
  const string filename = fullPath("rast100.sd");
  ifstream infile(filename.c_str(), ios::in | ios::binary);
  SurfData sd(infile, true);
  sd.write(fullPath("frombinary.txt"));
  infile.close();

  // Read it in as text as well and make sure the values are the same
  const string filenameText = fullPath("rast100.txt");
  ifstream infileText(filenameText.c_str(), ios::in);
  SurfData sdText(infileText, false);
  infileText.close();
  sdText.write(fullPath("fromtext.txt"));

  CPPUNIT_ASSERT_EQUAL(sd.points.size(), pointsInFile);
  CPPUNIT_ASSERT_EQUAL(sd.xsize, static_cast<unsigned>(2));
  CPPUNIT_ASSERT_EQUAL(sd.fsize, static_cast<unsigned>(1));
  CPPUNIT_ASSERT_EQUAL(sd.mapping.size(), pointsInFile);
  for (unsigned i = 0; i < pointsInFile; i++) {
    CPPUNIT_ASSERT_EQUAL(sd.mapping[i], i);
  }
  CPPUNIT_ASSERT(sd.excludedPoints.empty());
  CPPUNIT_ASSERT_EQUAL(sd.defaultIndex, unsignedZero);
  CPPUNIT_ASSERT_EQUAL(sd.xMatrix, static_cast<double*>(0));
  CPPUNIT_ASSERT_EQUAL(sd.yVector, static_cast<double*>(0));
  CPPUNIT_ASSERT(!sd.valid.xMatrix);
  CPPUNIT_ASSERT(!sd.valid.yVector);

  for (unsigned j = 0; j < pointsInFile; j++) {
    CPPUNIT_ASSERT_EQUAL(sd[j], sdText[j]);
  }
}

void SurfDataTest::testCopyConstructorSimple()
{
  unsigned pointsInFile = 100;
  const string filename = fullPath("rast100.sd");
  ifstream infile(filename.c_str(), ios::in);
  SurfData sd2(infile, true);
  infile.close();

  SurfData sd(sd2);
  
  CPPUNIT_ASSERT_EQUAL(sd.points.size(), pointsInFile);
  CPPUNIT_ASSERT_EQUAL(sd.xsize, static_cast<unsigned>(2));
  CPPUNIT_ASSERT_EQUAL(sd.fsize, static_cast<unsigned>(1));
  CPPUNIT_ASSERT_EQUAL(sd.mapping.size(), pointsInFile);
  for (unsigned i = 0; i < pointsInFile; i++) {
    CPPUNIT_ASSERT_EQUAL(sd.mapping[i], i);
  }
  CPPUNIT_ASSERT(sd.excludedPoints.empty());
  CPPUNIT_ASSERT_EQUAL(sd.defaultIndex, unsignedZero);
  CPPUNIT_ASSERT_EQUAL(sd.xMatrix, static_cast<double*>(0));
  CPPUNIT_ASSERT_EQUAL(sd.yVector, static_cast<double*>(0));
  CPPUNIT_ASSERT(!sd.valid.xMatrix);
  CPPUNIT_ASSERT(!sd.valid.yVector);

  // Compare the points, one by one
  // Assumes that SurfPoint::operator== works
  for (unsigned j = 0; j < pointsInFile; j++) {
    CPPUNIT_ASSERT_EQUAL(sd[j], sd2[j]);
  }
  
}

void SurfDataTest::testCopyConstructorComplex()
{
  unsigned pointsInFile = 100;
  const string filename = fullPath("rast100.sd");
  ifstream infile(filename.c_str(), ios::in);
  SurfData sd2(infile, true);
  infile.close();

  // Make some changes to object before copying
  // Change mapping
  skipPoints.insert(5); 
  skipPoints.insert(10); 
  skipPoints.insert(15); 
  unsigned numSkippedPoints = 3;
  sd2.setExcludedPoints(skipPoints);

  // Request XMatrix and YVector (results should be cached, which
  // will in turn change the results of the copy)
  const double* X = sd2.getXMatrix();
  const double* Y = sd2.getYVector();

  SurfData sd(sd2);

  
  
  CPPUNIT_ASSERT_EQUAL(sd.points.size(), pointsInFile);
  CPPUNIT_ASSERT_EQUAL(sd.mapping.size(), pointsInFile - numSkippedPoints);
  CPPUNIT_ASSERT_EQUAL(sd.size(), pointsInFile - numSkippedPoints);
  CPPUNIT_ASSERT_EQUAL(sd.xsize, static_cast<unsigned>(2));
  CPPUNIT_ASSERT_EQUAL(sd.fsize, static_cast<unsigned>(1));
  for (unsigned i = 0; i < sd.size(); i++) {
    CPPUNIT_ASSERT_EQUAL(sd.mapping[i], sd2.mapping[i]);
  }
  CPPUNIT_ASSERT_EQUAL(sd.excludedPoints.size(), 
    static_cast<unsigned>(numSkippedPoints));
  CPPUNIT_ASSERT_EQUAL(sd.defaultIndex, unsignedZero);
  CPPUNIT_ASSERT(sd.xMatrix != static_cast<double*>(0));
  CPPUNIT_ASSERT(sd.yVector != static_cast<double*>(0));
  CPPUNIT_ASSERT(sd.valid.xMatrix);
  CPPUNIT_ASSERT(sd.valid.yVector);

  // Compare the points, one by one
  // Assumes that SurfPoint::operator== works
  for (unsigned j = 0; j < sd2.size(); j++) {
    CPPUNIT_ASSERT_EQUAL(sd[j], sd2[j]);
  }
  
}

void SurfDataTest::testCopyActive()
{
  skipPoints.insert(1); 
  skipPoints.insert(4); 
  unsigned numSkippedPoints = 2;
  sdPtr1->setExcludedPoints(skipPoints);
  sdPtr1->setDefaultIndex(2);

  SurfData sd = sdPtr1->copyActive();
  CPPUNIT_ASSERT_EQUAL(sd.points.size(), numPoints - numSkippedPoints);
  CPPUNIT_ASSERT_EQUAL(sd.mapping.size(), numPoints - numSkippedPoints);
  CPPUNIT_ASSERT_EQUAL(sd.size(), numPoints - numSkippedPoints);
  CPPUNIT_ASSERT_EQUAL(sd.xsize, static_cast<unsigned>(3));
  CPPUNIT_ASSERT_EQUAL(sd.fsize, static_cast<unsigned>(3));
  CPPUNIT_ASSERT_EQUAL(sd.excludedPoints.size(), unsignedZero);
  CPPUNIT_ASSERT_EQUAL(sd.defaultIndex, static_cast<unsigned>(2));
  CPPUNIT_ASSERT(!sd.xMatrix );
  CPPUNIT_ASSERT(!sd.yVector );
  CPPUNIT_ASSERT(!sd.valid.xMatrix);
  CPPUNIT_ASSERT(!sd.valid.yVector);

  // Compare the points, one by one
  // Assumes that SurfPoint::operator== works
  for (unsigned j = 0; j < sd.size(); j++) {
    CPPUNIT_ASSERT_EQUAL(sd[j], (*sdPtr1)[j]);
  }
} 

void SurfDataTest::testCopyActiveEmpty()
{
  unsigned numSkippedPoints = skipAllPoints.size();
  sdPtr1->setExcludedPoints(skipAllPoints);

  SurfData sd = sdPtr1->copyActive();
  CPPUNIT_ASSERT_EQUAL(sd.points.size(), numPoints - numSkippedPoints);
  CPPUNIT_ASSERT_EQUAL(sd.mapping.size(), numPoints - numSkippedPoints);
  CPPUNIT_ASSERT_EQUAL(sd.size(), numPoints - numSkippedPoints);
  CPPUNIT_ASSERT_EQUAL(sd.xsize, static_cast<unsigned>(0));
  CPPUNIT_ASSERT_EQUAL(sd.fsize, static_cast<unsigned>(0));
  CPPUNIT_ASSERT_EQUAL(sd.excludedPoints.size(), unsignedZero);
  CPPUNIT_ASSERT_EQUAL(sd.defaultIndex, static_cast<unsigned>(0));
  CPPUNIT_ASSERT(!sd.xMatrix );
  CPPUNIT_ASSERT(!sd.yVector );
  CPPUNIT_ASSERT(!sd.valid.xMatrix);
  CPPUNIT_ASSERT(!sd.valid.yVector);
} 

void SurfDataTest::testAssignment()
{
  unsigned pointsInFile = 100;
  const string filename = fullPath("rast100.txt");
  ifstream infile(filename.c_str(), ios::in);
  SurfData sdText(infile, false);
  infile.close();

  SurfData sdMain(sdText);
  
  const string filenameBinary = fullPath("manypts.sd");
  ifstream infileBinary(filenameBinary.c_str(), ios::in);
  SurfData sdBinary(infileBinary, true);
  infileBinary.close();

  // Call assignment operator
  sdMain = sdBinary;
  
  CPPUNIT_ASSERT_EQUAL(sdMain.points.size(), sdBinary.size());
  CPPUNIT_ASSERT_EQUAL(sdMain.xsize, sdBinary.xSize());
  CPPUNIT_ASSERT_EQUAL(sdMain.fsize, sdBinary.fSize());
  CPPUNIT_ASSERT_EQUAL(sdMain.mapping.size(), sdBinary.mapping.size());
  for (unsigned i = 0; i < pointsInFile; i++) {
    CPPUNIT_ASSERT_EQUAL(sdMain.mapping[i], i);
  }
  CPPUNIT_ASSERT(sdMain.excludedPoints.empty());
  CPPUNIT_ASSERT_EQUAL(sdMain.defaultIndex, unsignedZero);
  CPPUNIT_ASSERT_EQUAL(sdMain.xMatrix, static_cast<double*>(0));
  CPPUNIT_ASSERT_EQUAL(sdMain.yVector, static_cast<double*>(0));
  CPPUNIT_ASSERT(!sdMain.valid.xMatrix);
  CPPUNIT_ASSERT(!sdMain.valid.yVector);

  // Compare the points, one by one
  // Assumes that SurfPoint::operator== works
  for (unsigned j = 0; j < sdBinary.size(); j++) {
    CPPUNIT_ASSERT_EQUAL(sdMain[j], sdBinary[j]);
  }
}

void SurfDataTest::testAssignmentToSelf()
{
  unsigned pointsInFile = 100;
  const string filename = fullPath("rast100.txt");
  ifstream infile(filename.c_str(), ios::in);
  SurfData sdText(infile, false);
  infile.close();

  SurfData sdMain(sdText);
  

  // Call assignment operator on self
  SurfData& sdRef = sdMain;
  sdMain = sdRef;
  
  CPPUNIT_ASSERT_EQUAL(sdMain.points.size(), pointsInFile);
  CPPUNIT_ASSERT_EQUAL(sdMain.xsize, static_cast<unsigned>(2));
  CPPUNIT_ASSERT_EQUAL(sdMain.fsize, static_cast<unsigned>(1));
  CPPUNIT_ASSERT_EQUAL(sdMain.mapping.size(), pointsInFile);
  for (unsigned i = 0; i < pointsInFile; i++) {
    CPPUNIT_ASSERT_EQUAL(sdMain.mapping[i], i);
  }
  CPPUNIT_ASSERT(sdMain.excludedPoints.empty());
  CPPUNIT_ASSERT_EQUAL(sdMain.defaultIndex, unsignedZero);
  CPPUNIT_ASSERT_EQUAL(sdMain.xMatrix, static_cast<double*>(0));
  CPPUNIT_ASSERT_EQUAL(sdMain.yVector, static_cast<double*>(0));
  CPPUNIT_ASSERT(!sdMain.valid.xMatrix);
  CPPUNIT_ASSERT(!sdMain.valid.yVector);
}

void SurfDataTest::testOperatorEquality()
{
  SurfData sd(*sdPtr1);
  CPPUNIT_ASSERT(sd == *sdPtr1);
  CPPUNIT_ASSERT(sd.operator==(*sdPtr1));
  CPPUNIT_ASSERT(sdPtr1->operator==(sd));

  SurfData sdText(fullPath("rast100.txt").c_str());
  SurfData sdBinary(fullPath("rast100.sd").c_str());
  CPPUNIT_ASSERT_EQUAL(sdText, sdBinary);
  sdText.points[0]->f[0] = 0.0;
  sdBinary.points[0]->f[0] = 1.0;
  CPPUNIT_ASSERT(!(sdText == sdBinary));
}

void SurfDataTest::testOperatorInequality()
{
  SurfData sd(*sdPtr1);
  CPPUNIT_ASSERT(!(sd != *sdPtr1));
  CPPUNIT_ASSERT(!sd.operator!=(*sdPtr1));
  CPPUNIT_ASSERT(!sdPtr1->operator!=(sd));

  SurfData sdText(fullPath("rast100.txt").c_str());
  SurfData sdBinary(fullPath("rast100.sd").c_str());
  CPPUNIT_ASSERT(!(sdText != sdBinary));
  sdText.points[0]->f[0] = 0.0;
  sdBinary.points[0]->f[0] = 1.0;
  CPPUNIT_ASSERT(sdText != sdBinary);
}

void SurfDataTest::testOperatorIndexing()
{
  const SurfPoint sp = (*sdPtr1)[0];
  CPPUNIT_ASSERT_EQUAL(sp, *sdPtr1->points[sdPtr1->mapping[0]]);
  skipPoints.insert(0);
  sdPtr1->setExcludedPoints(skipPoints);
  const SurfPoint sp2 = (*sdPtr1)[0];
  CPPUNIT_ASSERT_EQUAL(sp2, *sdPtr1->points[sdPtr1->mapping[0]]);
  CPPUNIT_ASSERT_EQUAL(sp2, *sdPtr1->points[1]);
}

void SurfDataTest::testOperatorIndexingBadIndex()
{
  // Should throw an exception because the index
  // is out of the range
  const SurfPoint sp = (*sdPtr1)[numPoints];
}

void SurfDataTest::testOperatorIndexingAnotherBadIndex()
{
  // Should throw an exception because the index
  // is out of the range
  const SurfPoint sp = (*sdPtr1)[static_cast<unsigned>(-1)];
}

void SurfDataTest::testSize()
{
  CPPUNIT_ASSERT_EQUAL(sdPtr1->size(), numPoints);
  skipPoints.insert(0);
  sdPtr1->setExcludedPoints(skipPoints);
  CPPUNIT_ASSERT_EQUAL(sdPtr1->size(), numPoints - 1);
  sdPtr1->setExcludedPoints(skipAllPoints);
  CPPUNIT_ASSERT_EQUAL(sdPtr1->size(), unsignedZero);
}

void SurfDataTest::testEmpty()
{
  CPPUNIT_ASSERT(!sdPtr1->empty());
  sdPtr1->setExcludedPoints(skipAllPoints);
  CPPUNIT_ASSERT(sdPtr1->empty());
}

void SurfDataTest::testXSize()
{
  CPPUNIT_ASSERT_EQUAL(sdPtr1->xSize(), dimPoints);
  surfpoints.clear();
  SurfData sd(surfpoints);
  CPPUNIT_ASSERT_EQUAL(sd.xSize(), unsignedZero);
  sdPtr1->setExcludedPoints(skipAllPoints);
  CPPUNIT_ASSERT_EQUAL(sdPtr1->xSize(), dimPoints);
}

void SurfDataTest::testFSize()
{
  CPPUNIT_ASSERT_EQUAL(sdPtr1->fSize(), numResponses);
  surfpoints.clear();
  SurfData sd(surfpoints);
  CPPUNIT_ASSERT_EQUAL(sd.fSize(), unsignedZero);
  sdPtr1->setExcludedPoints(skipAllPoints);
  CPPUNIT_ASSERT_EQUAL(sdPtr1->fSize(), numResponses);
}

void SurfDataTest::testGetExcludedPoints()
{
  const set<unsigned>& exPoints = sdPtr1->getExcludedPoints();
  CPPUNIT_ASSERT(exPoints.empty());
  sdPtr1->setExcludedPoints(skipAllPoints);
  const set<unsigned>& moreExPoints = sdPtr1->getExcludedPoints();
  CPPUNIT_ASSERT(moreExPoints == sdPtr1->excludedPoints);
}

void SurfDataTest::testGetResponse()
{
  CPPUNIT_ASSERT_EQUAL(sdPtr1->getResponse(0), surfpoints[0].f[0]);
  CPPUNIT_ASSERT_EQUAL(sdPtr1->getResponse(1), surfpoints[1].f[0]);
  CPPUNIT_ASSERT_EQUAL(sdPtr1->getResponse(2), surfpoints[2].f[0]);
  sdPtr1->setDefaultIndex(1);
  CPPUNIT_ASSERT_EQUAL(sdPtr1->getResponse(0), surfpoints[0].f[1]);
  CPPUNIT_ASSERT_EQUAL(sdPtr1->getResponse(1), surfpoints[1].f[1]);
  CPPUNIT_ASSERT_EQUAL(sdPtr1->getResponse(2), surfpoints[2].f[1]);
  sdPtr1->setDefaultIndex(2);
  CPPUNIT_ASSERT_EQUAL(sdPtr1->getResponse(2), surfpoints[2].f[2]);
  CPPUNIT_ASSERT_EQUAL(sdPtr1->getResponse(3), surfpoints[3].f[2]);
  CPPUNIT_ASSERT_EQUAL(sdPtr1->getResponse(4), surfpoints[4].f[2]);
  skipPoints.insert(3);
  sdPtr1->setExcludedPoints(skipPoints);
  CPPUNIT_ASSERT_EQUAL(sdPtr1->getResponse(1), surfpoints[1].f[2]);
  CPPUNIT_ASSERT_EQUAL(sdPtr1->getResponse(2), surfpoints[2].f[2]);
  CPPUNIT_ASSERT_EQUAL(sdPtr1->getResponse(3), surfpoints[4].f[2]);

}

void SurfDataTest::testGetDefaultIndex()
{
  CPPUNIT_ASSERT_EQUAL(sdPtr1->getDefaultIndex(), unsignedZero);
  sdPtr1->setDefaultIndex(1);
  CPPUNIT_ASSERT_EQUAL(sdPtr1->getDefaultIndex(), static_cast<unsigned>(1));
  sdPtr1->setExcludedPoints(skipAllPoints);
  CPPUNIT_ASSERT_EQUAL(sdPtr1->getDefaultIndex(), static_cast<unsigned>(1));
}

void SurfDataTest::testGetXMatrix()
{
  CPPUNIT_ASSERT(!sdPtr1->valid.xMatrix);
  const double* x = sdPtr1->getXMatrix();
  for (unsigned i = 0; i < numPoints; i++) {
    for (unsigned j = 0; j < dimPoints; j++) {
      CPPUNIT_ASSERT_EQUAL(x[i+j*numPoints], 
        sdPtr1->points[sdPtr1->mapping[i]]->x[j]);
    }
  }
  CPPUNIT_ASSERT(sdPtr1->valid.xMatrix);
  // Repeat function call after xMatrix has already been validated
  x = sdPtr1->getXMatrix();
  for (unsigned i = 0; i < numPoints; i++) {
    for (unsigned j = 0; j < dimPoints; j++) {
      CPPUNIT_ASSERT_EQUAL(x[i+j*numPoints], 
        sdPtr1->points[sdPtr1->mapping[i]]->x[j]);
    }
  }
  CPPUNIT_ASSERT(sdPtr1->valid.xMatrix);
}

void SurfDataTest::testGetYVector()
{
  CPPUNIT_ASSERT(!sdPtr1->valid.yVector);
  const double* y = sdPtr1->getYVector();
  for (unsigned i = 0; i < numPoints; i++) {
    CPPUNIT_ASSERT_EQUAL(y[i], 
      sdPtr1->points[sdPtr1->mapping[i]]->f[0]);
  }
  CPPUNIT_ASSERT(sdPtr1->valid.yVector);
  sdPtr1->setDefaultIndex(1);
  CPPUNIT_ASSERT(!sdPtr1->valid.yVector);
  y = sdPtr1->getYVector();
  CPPUNIT_ASSERT(sdPtr1->valid.yVector);
  for (unsigned i = 0; i < numPoints; i++) {
    CPPUNIT_ASSERT_EQUAL(y[i], 
      sdPtr1->points[sdPtr1->mapping[i]]->f[1]);
  }
  // repeat function call after yVector has been validated
  y = sdPtr1->getYVector();
  CPPUNIT_ASSERT(sdPtr1->valid.yVector);
  for (unsigned i = 0; i < numPoints; i++) {
    CPPUNIT_ASSERT_EQUAL(y[i], 
      sdPtr1->points[sdPtr1->mapping[i]]->f[1]);
  }
  CPPUNIT_ASSERT(sdPtr1->valid.yVector);
}

void SurfDataTest::testHasBinaryExtension()
{
  CPPUNIT_ASSERT(surfpack::hasExtension("abc.sd",".sd"));
  CPPUNIT_ASSERT(surfpack::hasExtension(".sd",".sd"));
  CPPUNIT_ASSERT(!surfpack::hasExtension(".sdx",".sd"));
  CPPUNIT_ASSERT(!surfpack::hasExtension("file.sdx",".sd"));
  CPPUNIT_ASSERT(!surfpack::hasExtension(".sd.txt",".sd"));
  CPPUNIT_ASSERT(surfpack::hasExtension("file.txt.sd",".sd"));
  CPPUNIT_ASSERT(surfpack::hasExtension("/dir/dir/dir/file.txt.sd",".sd"));
}

void SurfDataTest::testHasTextExtension()
{
  CPPUNIT_ASSERT(surfpack::hasExtension("abc.txt",".txt"));
  CPPUNIT_ASSERT(surfpack::hasExtension(".txt",".txt"));
  CPPUNIT_ASSERT(!surfpack::hasExtension(".txtx",".txt"));
  CPPUNIT_ASSERT(!surfpack::hasExtension(".tmt",".txt"));
  CPPUNIT_ASSERT(!surfpack::hasExtension("file.txx",".txt"));
  CPPUNIT_ASSERT(!surfpack::hasExtension(".txt.sd",".txt"));
  CPPUNIT_ASSERT(!surfpack::hasExtension("file. txt",".txt"));
  CPPUNIT_ASSERT(surfpack::hasExtension("file.sd.txt",".txt"));
  CPPUNIT_ASSERT(surfpack::hasExtension("/dir/dir/dir/file.sd.txt",".txt"));
}

// Commands

void SurfDataTest::testSetDefaultIndex()
{
  sdPtr1->setDefaultIndex(1);
  CPPUNIT_ASSERT(sdPtr1->defaultIndex == 1);
}

void SurfDataTest::testSetDefaultIndexBadIndex()
{
  sdPtr1->setDefaultIndex(numResponses);
}

void SurfDataTest::testSetDefaultIndexNoResponses()
{
  SurfPoint sp = surfpoints[0];
  sp.f.clear();
  vector<SurfPoint> pts;
  pts.push_back(sp);
  SurfData sd(pts);
  //Should throw an exception
  sd.setDefaultIndex(1);
}

void SurfDataTest::testSetResponse()
{
  CPPUNIT_ASSERT(sdPtr1->points[1]->f[0] != 0.0);
  sdPtr1->setResponse(1, -1.0);
  CPPUNIT_ASSERT(sdPtr1->points[1]->f[0] == -1.0);
  skipPoints.insert(1);
  skipPoints.insert(2);
  skipPoints.insert(3);
  sdPtr1->setDefaultIndex(1);
  sdPtr1->setExcludedPoints(skipPoints);
  sdPtr1->setResponse(0, -2.0);
  sdPtr1->setResponse(1, -3.0);
  skipPoints.clear();
  sdPtr1->setExcludedPoints(skipPoints);
  CPPUNIT_ASSERT(sdPtr1->getResponse(0) == -2.0);
  CPPUNIT_ASSERT(sdPtr1->getResponse(numPoints-1) == -3.0);
  sdPtr1->setDefaultIndex(0);
  CPPUNIT_ASSERT(sdPtr1->getResponse(1) == -1.0);
}

void SurfDataTest::testSetResponseBadIndex()
{
  // Should throw an exception
  sdPtr1->setResponse(numPoints, 0.0);
}

void SurfDataTest::testSetResponseAnotherBadIndex()
{
  sdPtr1->setExcludedPoints(skipAllPoints);
  // Should throw an exception
  sdPtr1->setResponse(numPoints, 0.0);
}

void SurfDataTest::testAddPoint()
{
  vector<double> x(dimPoints);
  vector<double> f(numResponses);
  SurfPoint sp(x,f);
  sdPtr1->addPoint(sp);
  CPPUNIT_ASSERT(sdPtr1->size() == numPoints+1);
  CPPUNIT_ASSERT((*sdPtr1)[numPoints] == sp);
  // skipAllPoints holds indices 0-4 but not 5 (yet)
  sdPtr1->setExcludedPoints(skipAllPoints);
  CPPUNIT_ASSERT(sdPtr1->size() == 1);
  x[0] = 1.0;
  SurfPoint sp2(x,f);
  sdPtr1->addPoint(sp2);
  CPPUNIT_ASSERT(sdPtr1->size() == 2);
  skipAllPoints.insert(numPoints);
  skipAllPoints.insert(numPoints+1);
  sdPtr1->setExcludedPoints(skipAllPoints);
  CPPUNIT_ASSERT(sdPtr1->empty());
  x[0] = 2.0;
  SurfPoint sp3(x,f);
  sdPtr1->addPoint(sp3);
  CPPUNIT_ASSERT(sdPtr1->size() == 1);
  sdPtr1->setExcludedPoints(skipPoints);
  CPPUNIT_ASSERT((*sdPtr1)[numPoints+1] == sp2);
  CPPUNIT_ASSERT((*sdPtr1)[numPoints+2] == sp3);
  CPPUNIT_ASSERT(!sdPtr1->valid.xMatrix);
  CPPUNIT_ASSERT(!sdPtr1->valid.yVector);
}

void SurfDataTest::testAddPointBadDimension()
{
  vector<double> x(dimPoints+1);
  vector<double> f(numResponses);
  SurfPoint sp(x,f);
  sdPtr1->addPoint(sp);
}

void SurfDataTest::testAddPointBadNumResponses()
{
  vector<double> x(dimPoints);
  vector<double> f(numResponses-1);
  SurfPoint sp(x,f);
  sdPtr1->addPoint(sp);
}

void SurfDataTest::testAddPointToEmptySet()
{
  vector<SurfPoint> emptySet;
  SurfData sd(emptySet);
  CPPUNIT_ASSERT(sd.xSize() == 0);
  CPPUNIT_ASSERT(sd.fSize() == 0);
  CPPUNIT_ASSERT(sd.empty());
  const double* x = sd.getXMatrix();
  const double* y = sd.getYVector();
  sd.addPoint(surfpoints[0]);
  CPPUNIT_ASSERT(sd.xSize() == dimPoints);
  CPPUNIT_ASSERT(sd.fSize() == numResponses);
  CPPUNIT_ASSERT(sd.size() == 1);
  CPPUNIT_ASSERT(!sd.valid.xMatrix);
  CPPUNIT_ASSERT(!sd.valid.yVector);
}

void SurfDataTest::testAddResponse()
{
  vector<double> newValues(numPoints, -1.0);
  sdPtr1->addResponse(newValues);
  CPPUNIT_ASSERT(sdPtr1->xSize() == dimPoints); 
  CPPUNIT_ASSERT(sdPtr1->fSize() == numResponses + 1); 
  sdPtr1->setDefaultIndex(numResponses);
  for (unsigned i = 0; i < numPoints; i++) {
    CPPUNIT_ASSERT(sdPtr1->getResponse(i) == -1.0);
  }
}

void SurfDataTest::testAddResponseToEmptySet()
{
  vector<SurfPoint> emptySet;
  SurfData sd(emptySet);
  vector<double> newValues(numPoints, -1.0);
  sd.addResponse(newValues);
}

void SurfDataTest::testAddResponseWithSkipped()
{
  sdPtr1->setExcludedPoints(skipAllPoints);
  vector<double> newValues(numPoints, -1.0);
  sdPtr1->addResponse(newValues);
}

void SurfDataTest::testAddResponseWrongNumber()
{
  vector<double> newValues(numPoints - 1, -1.0);
  sdPtr1->addResponse(newValues);
}

void SurfDataTest::testSetExcludedPoints()
{
  skipPoints.insert(1);
  skipPoints.insert(3);
  sdPtr1->setExcludedPoints(skipPoints);
  CPPUNIT_ASSERT(sdPtr1->excludedPoints == skipPoints);
  CPPUNIT_ASSERT(sdPtr1->mapping.size() == 3);
  CPPUNIT_ASSERT(sdPtr1->mapping[0] == 0);
  CPPUNIT_ASSERT(sdPtr1->mapping[1] == 2);
  CPPUNIT_ASSERT(sdPtr1->mapping[2] == 4);
  skipPoints.clear();
  skipPoints.insert(0);
  skipPoints.insert(4);
  sdPtr1->setExcludedPoints(skipPoints);
  CPPUNIT_ASSERT(sdPtr1->excludedPoints == skipPoints);
  CPPUNIT_ASSERT(sdPtr1->mapping.size() == 3);
  CPPUNIT_ASSERT(sdPtr1->mapping[0] == 1);
  CPPUNIT_ASSERT(sdPtr1->mapping[1] == 2);
  CPPUNIT_ASSERT(sdPtr1->mapping[2] == 3);
}

void SurfDataTest::testSetExcludedPointsToNone()
{
  sdPtr1->setExcludedPoints(skipAllPoints);
  sdPtr1->setExcludedPoints(skipPoints);
  CPPUNIT_ASSERT(sdPtr1->excludedPoints.empty());
  CPPUNIT_ASSERT(sdPtr1->mapping.size() == numPoints);
}

void SurfDataTest::testSetExcludedPointsTooMany()
{
  skipAllPoints.insert(numPoints+1);
  sdPtr1->setExcludedPoints(skipAllPoints);
}

void SurfDataTest::testWriteBinary()
{
  sdPtr1->write(fullPath("surfdata1.sd"));
  SurfData sd(fullPath("surfdata1.sd"));
  CPPUNIT_ASSERT(*sdPtr1 == sd);
}

void SurfDataTest::testWriteText()
{
  sdPtr1->write(fullPath("surfdata1.txt"));
  SurfData sd(fullPath("surfdata1.txt"));
  CPPUNIT_ASSERT(*sdPtr1 == sd);
}

void SurfDataTest::testWriteNoPoints()
{
  sdPtr1->setExcludedPoints(skipAllPoints);
  sdPtr1->write(fullPath("surfdata1.txt"));
}

void SurfDataTest::testWriteNoFile()
{
  // This may be a platform-dependent test.
  // A better test might be to try to open a file
  // when the user does not have write privileges
  sdPtr1->write("///.txt");
}

void SurfDataTest::testWriteBadFileExtension()
{
  sdPtr1->write(fullPath("file.zip"));
}

void SurfDataTest::testReadNoFile()
{
  sdPtr1->read(fullPath("file_does_not_exist.txt"));
}

void SurfDataTest::testReadTextFileTooShort()
{
  sdPtr1->read(fullPath("claimsTooMany.txt"));
}

void SurfDataTest::testReadBinaryFileTooShort()
{
  sdPtr1->read(fullPath("claimsTooMany.sd"));
}

void SurfDataTest::testStreamInsertion()
{
  // If this test doesn't throw an exception,
  // it is presumed to have worked
  ofstream blackhole("/dev/null",ios::out);
  blackhole << (*sdPtr1) << endl;
  cout << "End SurfDataTest" << endl;
}

