// These lines were automatically prepended
#ifndef SURFPACK_CONFIG_H
#define SURFPACK_CONFIG_H
#include "config.h"
#endif
// End of prepended lines

// ----------------------------------------------------------
// Project: SURFPACK++
//
// File:        SurfData.cpp
// Author:      Eric Cyr
// Modified:	Mark Richards
//
// Description: 
// + SurfData class - this is the container for all the data 
//   from which an empirical model is created. 
// + Left shift (<<) operator for SurfData. 
// ----------------------------------------------------------

#include <vector>
#include <set>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>

#include "SurfData.h"
#include "SurfPoint.h"
//#include "Surface.h"

using namespace std;

// ____________________________________________________________________________
// Creation, Destruction, Initialization 
// ____________________________________________________________________________

/// Vector of points will be copied
SurfData::SurfData(const vector<SurfPoint>& points) : valid()
{
#ifdef __TESTING_MODE__
  constructCount++;
#endif
  if (points.empty()) {
    this->xsize = 0;
    this->fsize = 0;
  } else {
    this->xsize = points[0].xSize();
    this->fsize = points[0].fSize();
    this->points = points;
  }
  init();
  // Check to make sure data points all have the same number of dimensions
  // and response values.  An exception will be thrown otherwise.
  sanityCheck();
}

/// Read a set of SurfPoints from a file
SurfData::SurfData(const string filename) : valid() 
{
#ifdef __TESTING_MODE__
  constructCount++;
#endif
  read(filename);
  init();
}
  
/// Read a set of SurfPoints from an istream
SurfData::SurfData(std::istream& is, bool binary) : valid()
{
  if (binary) {
    readBinary(is);
  } else {
    readText(is);
  }
  init();
}


/// Makes a deep copy of the object 
SurfData::SurfData(const SurfData& other) : xsize(other.xsize), 
  fsize(other.fsize), points(other.points),excludedPoints(other.excludedPoints),
  mapping(other.mapping), valid(other.valid)
{
#ifdef __TESTING_MODE__
  constructCount++;
  copyCount++;
#endif
  copyBlockData(other);
}

/// STL data members' resources automatically deallocated 
SurfData::~SurfData() 
{
  // clean up the related surfaces
  //listeners.clear();
  // clean up all the created points
  //deletePoints();
#ifdef __TESTING_MODE__
  destructCount++;
#endif
  cleanup();
}

void SurfData::init()
{
  defaultIndex = 0;
  defaultMapping();
  xMatrix = 0;
  yVector = 0;
}

/// Copy only the "active" points
SurfData* SurfData::copyActive()
{
  vector<SurfPoint> activePoints;
  for (unsigned i = 0; i < mapping.size(); i++) {
    activePoints.push_back(points[mapping[i]]);
  }
  SurfData* newSD = new SurfData(activePoints);
  if (!activePoints.empty()) {
    newSD->setDefaultIndex(defaultIndex);
  }
  return newSD;
}
  
// Copy xMatrix and yVector from another SurfData object
void SurfData::copyBlockData(const SurfData& other)
{
  try {
    if (other.valid.xMatrix) {
      unsigned numElements = mapping.size() * xsize;
      xMatrix = new double[numElements];
      memcpy(xMatrix,other.xMatrix, numElements*sizeof(double));
    } else {
      xMatrix = 0;
    }
 
    if (other.valid.yVector) {
      unsigned numElements = mapping.size();
      yVector = new double[numElements];
      memcpy(yVector,other.yVector, numElements*sizeof(double));
    } else {
      yVector = 0;
    }
  } catch (...) {
    cerr << "Error in SurfData::copyBlockData.  Rethrowing exception." << endl;
    // If xMatrix allocation succeeded but yVector failed, cleanup necessary
    cleanup();
    throw;
  }
}
  
// Deallocate any memory allocated for xMatrix and/or yVector
void SurfData::cleanup()
{
  delete xMatrix;
  delete yVector;
}
// ____________________________________________________________________________
// Overloaded operators 
// ____________________________________________________________________________

/// makes deep copy 
SurfData& SurfData::operator=(const SurfData& sd)
{
  if (*this != sd) {
    this->xsize = sd.xsize;
    this->fsize = sd.fsize;
    this->points = sd.points;
    this->excludedPoints = sd.excludedPoints;
    this->mapping = sd.mapping;
    this->valid = sd.valid;
    this->defaultIndex = sd.defaultIndex;
    cleanup();
    copyBlockData(sd);

    // now all the surfaces need to be update
    // or at least told of a change
    //notifyListeners();
  }
  return (*this);
}

/// makes deep comparison
bool SurfData::operator==(const SurfData& sd)
{
  return (this->xsize == sd.xsize && 
          this->fsize == sd.fsize &&
          this->size() == sd.size() &&
          this->points == sd.points);
}
      
/// makes deep comparison
bool SurfData::operator!=(const SurfData& sd)
{
  return !(*this == sd);
}

/// Return a point from the data set
//SurfPoint& SurfData::Point(unsigned index) 
SurfPoint& SurfData::operator[](unsigned index) 
{
  static string header("Indexing error in SurfData::operator[].");
  checkRange(header, index);
  return points[mapping[index]];
}

/// Return a point from the data set
//SurfPoint& SurfData::Point(unsigned index) 
const SurfPoint& SurfData::operator[](unsigned index) const
{
  static string header("Indexing error in SurfData::operator[] const.");
  checkRange(header, index);
  return points[mapping[index]];
}

// ____________________________________________________________________________
// Queries 
// ____________________________________________________________________________

/// Return the number of SurfPoints in the data set 
unsigned SurfData::size() const 
{ 
  return mapping.size(); 
}

/// True if there are no points
bool SurfData::empty() const
{
  return mapping.empty();
}

/// Return the dimensionality of the SurfPoints 
unsigned SurfData::xSize() const 
{ 
  return xsize; 
}

/// Return the number of response values in the data set
unsigned SurfData::fSize() const 
{ 
  return fsize; 
}

/// Return the set of excluded points (the indices)
const std::set<unsigned>& SurfData::getExcludedPoints() const 
{
  return excludedPoints;
}

/// Get the response value of the (index)th point that corresponds to this
/// surface
double SurfData::getResponse(unsigned index) const
{
  static string header("Indexing error in SurfData::getResponse.");
  checkRange(header, index);
  return points[mapping[index]].F(defaultIndex);
}

/// Get default index
unsigned SurfData::getDefaultIndex() const
{
  return defaultIndex;
}

/// Return a reference to the SurfPoints vector 
//vector<SurfPoint>& SurfData::Points() 
//{ 
//  return points; 
//}

/// Return point domains as a matrix in a contiguous block.  Be careful.
/// The data should not be changed.
const double* SurfData::getXMatrix() const
{
  if (!valid.xMatrix) {
    validateXMatrix();
  }
  return xMatrix;
}

/// Return response values for the default response in a contiguous block.  
/// Be careful. The data should not be changed.
const double* SurfData::getYVector() const
{
  if (!valid.yVector) {
    validateYVector();
  }
  return yVector;
}
  
/// Returns true if the filename extension is .sd.
bool SurfData::hasBinaryExtension(const std::string& filename)
{
  return (filename.find(".sd") == filename.size() - 3);
}

/// Returns true if the filename extension is .txt.
bool SurfData::hasTextExtension(const std::string& filename)
{
  return (filename.find(".txt") == filename.size() - 4);
}
  

// ____________________________________________________________________________
// Commands 
// ____________________________________________________________________________

/// Specify which response value getResponse will return
void SurfData::setDefaultIndex(unsigned index)
{
  static string header("Indexing error in SurfData::setDefaultIndex.");
  checkRange(header, index);
  valid.yVector = false;
  defaultIndex = index;
}
  
/// Set the response value of the (index)th point that corresponds to this
/// surface
void SurfData::setResponse(unsigned index, double value)
{
  static string header("Indexing error in SurfData::setResponse.");
  checkRange(header, index);
  points[mapping[index]].F(defaultIndex, value);
  valid.yVector = false;
}
  
/// Add a point to the data set. The parameter point will be copied.
void SurfData::addPoint(const SurfPoint& sp) 
{
  if (points.empty()) {
    xsize = sp.xSize();
    fsize = sp.fSize();
  } else {
    if (sp.xSize() != xsize || sp.fSize() != fsize) {
      ostringstream errormsg;
      errormsg << "Error in SurfData::addPoint.  Points in this data set "
	       << "have " << xsize << " dimensions and " << fsize
	       << " response values; point to be added has "
	       << sp.xSize() << " dimensions and " << sp.fSize()
	       << " response values." << endl;
      throw bad_surf_data(errormsg.str());
    }
  }
  points.push_back(sp);
  mapping.push_back(points.size()-1);
  valid.xMatrix = false;
  valid.yVector = false;
}

/// Add a new response variable to each point. 
/// Return the index of the new variable.
unsigned SurfData::addResponse(const vector<double>& newValues)
{
  unsigned newindex;
  ostringstream errormsg;
  if (points.empty()) {
    throw bad_surf_data(
             "Cannot add response because there are no data points"
          );
  } else if (points.size() != mapping.size()) {
    errormsg << "Cannot add response because physical set size is different "
	     << "than logical set size.\nBefore adding another response, "
             << "clear \"excluded points\" or create a new data set by using " 
	     << "the SurfData::copyActive method." << endl;
    throw bad_surf_data(errormsg.str());
  } else if (newValues.size() != points.size()) {
    errormsg << "Cannot add another response because the number of new response"
             << " values does not match the size of the physical data set." 
             << endl;
    throw bad_surf_data(errormsg.str());
  } else {
    newindex = points[mapping[0]].addResponse(newValues[0]);
    fsize++;
    for (unsigned i = 1; i < points.size(); i++) {
      newindex = points[mapping[i]].addResponse(newValues[i]);
      assert(newindex == fsize - 1);
    }
  }
  return newindex;
}

/// Specify which points should be skipped
void SurfData::setExcludedPoints(const std::set<unsigned>& excludedPoints)
{
  if (excludedPoints.size() > points.size()) {
    throw bad_surf_data(
      "Size of set of excluded points exceeds size of SurfPoint set"
    );
  } else if (excludedPoints.empty()) {
    defaultMapping();
    this->excludedPoints.clear();
  } else {
    // The size of the logical data set is the size of the physical
    // data set less the number of excluded points    
    mapping.resize(points.size() - excludedPoints.size());
    unsigned mappingIndex = 0;
    unsigned sdIndex = 0;
    // map the valid indices to the physical points in points
    while (sdIndex < points.size()) {
      if (excludedPoints.find(sdIndex) == excludedPoints.end()) {
        mapping[mappingIndex++] = sdIndex;
      }
      sdIndex++;
    }
    this->excludedPoints = excludedPoints;
    assert(mappingIndex == mapping.size());
  }
}

//void SurfData::addListener(Surface * surface)
//{
//  /// only add the listener if its not already there
//  list<Surface*>::iterator itr = 
//    find(listeners.begin(),listeners.end(),surface);
//  if(itr==listeners.end()) {
//    listeners.push_back(surface);
//  }
//}
//
//void SurfData::removeListener(Surface * surface)
//{
//  // make sure its OK to erase the object, then do so
//  list<Surface*>::iterator itr =
//    find(listeners.begin(),listeners.end(),surface);
//  if(itr!=listeners.end()) {
//    listeners.erase(itr);
//  }
//}

/// Maps all indices to themselves
void SurfData::defaultMapping()
{
  mapping.resize(points.size());
  for (unsigned i = 0; i < points.size(); i++) {
    mapping[i] = i;
  }
}
   
/// Creates a matrix of the domains for all of the points in a contiguous
/// block of memory, for use in matrix operations
void SurfData::validateXMatrix() const
{
  delete xMatrix;
  xMatrix = new double[mapping.size() * xsize];
  for (unsigned pt = 0; pt < mapping.size(); pt++) {
    for (unsigned xval = 0; xval < xsize; xval++) {
      xMatrix[xval+pt*mapping.size()] = points[mapping[pt]].X()[xval];
    }
  }
  valid.xMatrix = true;
}

/// Creates a vector of response values for the default response value in
/// a contiguous blcok of memory
void SurfData::validateYVector() const
{
  delete yVector;
  yVector = new double[mapping.size()];
  for (unsigned pt = 0; pt < mapping.size(); pt++) {
    yVector[pt] = getResponse(pt);
  }
  valid.yVector = true;
}

// ____________________________________________________________________________
// I/O
// ____________________________________________________________________________

/// Write a set of SurfPoints to a file.  Opens the file and calls other version
void SurfData::write(const std::string& filename) const
{
  if (mapping.empty()) {
    ostringstream errormsg;
    errormsg << "Cannot write SurfData object to stream."
	     << "  No active data points." << endl;
    throw bad_surf_data(errormsg.str());
  }
  bool binary = testFileExtension(filename);
  ofstream outfile(filename.c_str(), (binary ? ios::out|ios::binary : ios::out));
  if (!outfile) {
    throw file_open_failure(filename);
  } else if (binary) {
    writeBinary(outfile);
  } else {
    writeText(outfile);
  }
  outfile.close();
}

/// Read a set of SurfPoints from a file.  Opens file and calls other version.
void SurfData::read(const string& filename)
{
  bool binary = testFileExtension(filename);
  // Open file in binary or text mode based on filename extension (.sd or .txt)
  ifstream infile(filename.c_str(), (binary ? ios::in|ios::binary : ios::in));
  if (!infile) {
    throw file_open_failure(filename);
  } else if (binary) {
    readBinary(infile);
  } else {
    readText(infile);
  }
  infile.close();
}

/// Write a set of SurfPoints to an output stream
void SurfData::writeBinary(ostream& os) const 
{
  try {
    unsigned s = mapping.size();
    os.write((char*)&s,sizeof(s));
    os.write((char*)&xsize,sizeof(xsize));
    os.write((char*)&fsize,sizeof(fsize));
    for (unsigned i = 0; i < mapping.size(); i++) {
      points[mapping[i]].writeBinary(os);
    }
  } catch(...) {
    cerr << "Unknown exception caught in SurfData::writeBinary" << endl;
    throw;
  }
}

/// Write a set of SurfPoints to an output stream
void SurfData::writeText(ostream& os) const
{
  try {
    // Write an extra space after last token to prevent
    // ifstream from setting the bad_bit when data is later read in.
    os << mapping.size() << " " << endl
       << xsize << " " << endl 
       << fsize << " " << endl;
    for (unsigned i = 0; i < mapping.size(); i++) {
      points[mapping[i]].writeText(os);
    }
  } catch(...) {
    cerr << "Unknown exception caught in SurfData::writeText" << endl;
    throw;
  }
}

/// Read a set of SurfPoints from an input stream
void SurfData::readBinary(istream& is) 
{
  try {
    unsigned size;  
    is.read((char*)&size,sizeof(size));
    is.read((char*)&xsize,sizeof(xsize));
    is.read((char*)&fsize,sizeof(fsize));
    points.clear();
    for (unsigned i = 0; i < size; i++) {
      // True for second argument signals a binary read
      checkForEOF(is);
      points.push_back(SurfPoint(xsize,fsize,is,true));  
    }
    defaultMapping();
  } catch(...) {
    cerr << "Unknown exception caught in SurfData::readBinary" << endl;
    throw;
  }
}

/// Read a set of SurfPoints from an input stream
void SurfData::readText(istream& is) 
{
  try {
    string sline;
    getline(is,sline);
    istringstream streamline(sline);
    unsigned size;
    streamline >> size;
    getline(is,sline);
    streamline.str(sline);
    streamline >> xsize;
    getline(is,sline);
    streamline.str(sline);
    streamline >> fsize;
    points.clear();
    for (unsigned i = 0; i < size; i++) {
      checkForEOF(is);
      // False for last argument signals a text read
      points.push_back(SurfPoint(xsize,fsize,is,false));  
    }
    defaultMapping();
  } catch(...) {
    cerr << "Unknown exception caught in SurfData::readText" << endl;
    throw;
  }
}

// so a SurfData object can be printed
ostream& operator<<(ostream& os, const SurfData& sd) 
{ 
  sd.writeText(os); 
  return os;
}

// ____________________________________________________________________________
// Helper methods 
// ____________________________________________________________________________

/// Returns true if file is opened in binary mode, false if it is opened
/// in text mode.  If file cannot be opened, an exception is thrown.
bool SurfData::testFileExtension(const std::string& filename) const
{
  if (hasBinaryExtension(filename)) {
    return true;
  } else if (hasTextExtension(filename)) {
    return false;
  } else {
    throw io_exception(
      "Unrecognized filename extension in SurfData::read.  Use .sd or .txt"
    );
  }
}

//void SurfData::notifyListeners()
//{
//  list<Surface*>::iterator itr;
//  for(itr=listeners.begin();itr!=listeners.end();++itr)
//    (*itr)->notify();
//}

// ____________________________________________________________________________
// Testing 
// ____________________________________________________________________________

// Throw an exception if there are any mismatches in the number of
// dimensions or number of response values among points in the data set  
void SurfData::sanityCheck() const
{
  if (!points.empty()) {
    unsigned dimensionality = points[0].xSize();
    unsigned numResponses = points[0].fSize();
    for (unsigned i = 1; i < points.size(); i++) {
      if (points[i].xSize() != dimensionality ||
          points[i].fSize() != numResponses ) {
        ostringstream errormsg;
        errormsg << "Error in SurfData::sanityCheck." << endl
		 << "Point 0 has " << dimensionality << " dimensions "
                 << "and " << numResponses << " response values, " << endl
                 << "but point " << i << " has " << points[i].xSize()
 		 << " dimensions and " << points[i].fSize() << "response "
 		 << " values.";
	throw bad_surf_data(errormsg.str());
      } // if mismatch
    } // for each point
  } // if !empty
}

/// Check that the index falls within acceptable boundaries (i.e., is
/// less than mapping.size()
void SurfData::checkRange(const string& header, unsigned index) const
{
  if (index > mapping.size()) {
    ostringstream errormsg;
    errormsg << header << endl;
    if (mapping.empty()) {
      errormsg << "Index " << index << " specified, but there are zero points "
	       << "in the logical data set."
               << endl;
    } else {
      errormsg << "Requested: " 
	     << index 
	     << "; actual max index: "
	     << mapping.size() - 1
	     << endl;
    }
    throw range_error(errormsg.str());
  }
}

/// Make sure eof has not been reached unexpectedly
void SurfData::checkForEOF(istream& is) const
{
  if (is.eof()) {
    throw io_exception("End of file reached unexpectedly.");
  }
}

void SurfData::writeMatrix(const string header, double* mat, unsigned rows, 
  unsigned columns, ostream& os)
{
  if (header != "none" && header != "") {
    cout << header << endl;
  }
  for (unsigned r = 0; r < rows; r++) {
    for (unsigned c = 0; c < columns; c++) {
      os << setw(15) << mat[r + c*rows];
    }
    os << endl;
  }
}

void SurfData::writeMatrix(const string filename, double* mat, unsigned rows, 
  unsigned columns)
{
  ofstream outfile(filename.c_str(),ios::out);
  if (!outfile) {
    cerr << "Could not open file (" << filename << ") in writeMatrix." << endl;
    return;
  }
  writeMatrix("none",mat,rows,columns,outfile);
  outfile.close();
}
#ifdef __TESTING_MODE__
  int SurfData::constructCount = 0;
  int SurfData::copyCount = 0;
  int SurfData::destructCount = 0;
#endif

