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
#include <cassert>

#include "surfpack.h"
#include "SurfData.h"
#include "SurfPoint.h"
#include "Surface.h"

using namespace std;
// ____________________________________________________________________________
// Constants 
// ____________________________________________________________________________
const int SurfData::GOING_OUT_OF_EXISTENCE = 1;
const int SurfData::DATA_MODIFIED = 2;

// ____________________________________________________________________________
// Creation, Destruction, Initialization 
// ____________________________________________________________________________

/// Vector of points will be copied
SurfData::SurfData(const vector<SurfPoint>& points_) : valid()
{
#ifdef __TESTING_MODE__
  constructCount++;
#endif
  init();
  if (points_.empty()) {
    this->xsize = 0;
    this->fsize = 0;
  } else {
    this->xsize = points_[0].xSize();
    this->fsize = points_[0].fSize();
    for (unsigned i = 0; i < points_.size(); i++) {
      this->addPoint(points_[i]);
    }
  }
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
  init();
  read(filename);
}
  
/// Read a set of SurfPoints from an istream
SurfData::SurfData(std::istream& is, bool binary) : valid()
{
  init();
  if (binary) {
    readBinary(is);
  } else {
    readText(is);
  }
}


/// Makes a deep copy of the object 
SurfData::SurfData(const SurfData& other) : xsize(other.xsize), 
  fsize(other.fsize),excludedPoints(other.excludedPoints),
  defaultIndex(other.defaultIndex)
{
#ifdef __TESTING_MODE__
  constructCount++;
  copyCount++;
#endif
  for (unsigned i = 0; i < other.points.size(); i++) {
    this->addPoint(*other.points[i]);
  }
  copyBlockData(other);
  valid = other.valid;
  mapping = other.mapping;
}

/// STL data members' resources automatically deallocated 
SurfData::~SurfData() 
{
  cleanupStarted = true;
  //cout << "Starting to cleanup: " << this << endl;
  notifyListeners(GOING_OUT_OF_EXISTENCE);
  // clean up the related surfaces
  listeners.clear();
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
  cleanupStarted = false;
}

/// Copy only the "active" points
SurfData SurfData::copyActive()
{
  vector<SurfPoint> activePoints;
  for (unsigned i = 0; i < mapping.size(); i++) {
    activePoints.push_back(*points[mapping[i]]);
  }
  SurfData newSD(activePoints);
  if (!activePoints.empty()) {
    newSD.setDefaultIndex(defaultIndex);
  }
  return newSD;
}
  
// Copy xMatrix and yVector from another SurfData object
void SurfData::copyBlockData(const SurfData& other)
{
  // Prepared to catch an exception if memory allocation fails
  //try {
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
}
  
// Deallocate any memory allocated for xMatrix and/or yVector
void SurfData::cleanup()
{
  delete xMatrix;
  delete yVector;
  xMatrix=0;
  yVector=0;
  valid.xMatrix = false;
  valid.yVector = false;
  mapping.clear();
  orderedPoints.clear();
  for (unsigned j = 0; j < points.size(); j++) {
    delete points[j];
    points[j] =0;
  }
  points.clear();
  excludedPoints.clear();
}
// ____________________________________________________________________________
// Overloaded operators 
// ____________________________________________________________________________

/// makes deep copy 
SurfData& SurfData::operator=(const SurfData& other)
{
  if (*this != other) {
    cleanup();
    this->xsize = other.xsize;
    this->fsize = other.fsize;
    for (unsigned i = 0; i < other.points.size(); i++) {
      this->addPoint(*other.points[i]);
    }
    this->excludedPoints = other.excludedPoints;
    this->mapping = other.mapping;
    this->valid = other.valid;
    this->defaultIndex = other.defaultIndex;
    copyBlockData(other);

  }
  return (*this);
}

/// makes deep comparison
bool SurfData::operator==(const SurfData& other) const
{
  if (this->xsize == other.xsize && 
      this->fsize == other.fsize &&
      this->size() == other.size()) {
    for (unsigned i = 0; i < points.size(); i++) {
      if (*this->points[i] != *other.points[i]) {
        return false;
      }
    }
    return true;
  } else {
    return false;
  }
}
      
/// makes deep comparison
bool SurfData::operator!=(const SurfData& other) const
{
  return !(*this == other);
}

/// Return a point from the data set
//SurfPoint& SurfData::Point(unsigned index) 
const SurfPoint& SurfData::operator[](unsigned index) const
{
  static string header("Indexing error in SurfData::operator[] const.");
  checkRangeNumPoints(header, index);
  return *points[mapping[index]];
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
  checkRangeNumPoints(header, index);
  return points[mapping[index]]->F(defaultIndex);
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
//bool SurfData::hasBinaryExtension(const std::string& filename)
//{
//  return (filename.find(".sd") == filename.size() - 3);
//}
//
///// Returns true if the filename extension is .txt.
//bool SurfData::hasTextExtension(const std::string& filename)
//{
//  return (filename.find(".txt") == filename.size() - 4);
//}
  

// ____________________________________________________________________________
// Commands 
// ____________________________________________________________________________

/// Specify which response value getResponse will return
void SurfData::setDefaultIndex(unsigned index)
{
  static string header("Indexing error in SurfData::setDefaultIndex.");
  checkRangeNumResponses(header, index);
  valid.yVector = false;
  defaultIndex = index;
}
  
/// Set the response value of the (index)th point that corresponds to this
/// surface
void SurfData::setResponse(unsigned index, double value)
{
  static string header("Indexing error in SurfData::setResponse.");
  checkRangeNumPoints(header, index);
  points[mapping[index]]->F(defaultIndex, value);
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
  SurfPointSet::iterator iter;
  // This should be a safe const cast.  All that's happening is a check
  // to see if another data point at the same location in the space has already
  // been added.
  iter = orderedPoints.find(const_cast<SurfPoint*>(&sp));
  if (iter == orderedPoints.end()) {
    //cout << "Checked for dup" << endl;
    points.push_back(new SurfPoint(sp));
    int beforesize = orderedPoints.size();
    orderedPoints.insert(points[points.size()-1]);
    int aftersize = orderedPoints.size();
    if (beforesize == aftersize) { 
      cout << "Something went wrong" << endl;
      cout << "beforesize: " << beforesize << " aftersize: " << aftersize << endl;
    }
    if (orderedPoints.size() != points.size()) {
      cout << "oPsize: " << orderedPoints.size() << " points.size(): " << points.size() << endl;
    }
    mapping.push_back(points.size()-1);
  } else {
    cerr << "Duplication" << endl;
    // Replace the old point with this new one
    SurfPoint* spPtr = *iter;
    *spPtr = sp;
  }
  valid.xMatrix = false;
  valid.yVector = false;
  notifyListeners(SurfData::DATA_MODIFIED);
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
    newindex = points[mapping[0]]->addResponse(newValues[0]);
    fsize++;
    for (unsigned i = 1; i < points.size(); i++) {
      newindex = points[mapping[i]]->addResponse(newValues[i]);
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

void SurfData::addListener(Surface* surface)
{
  /// only add the listener if its not already there
  cout << "SurfData: " << this
       << " #listbef: " << listeners.size()
       << " adding: " << surface;
  list<Surface*>::iterator itr = 
    find(listeners.begin(),listeners.end(),surface);
  if(itr ==listeners.end() ) {
    listeners.push_back(surface);
    //cout << "Listener added: " << listeners.size() 
    // 	   << "this: " << this << endl;
  }
  cout << " #listaft: " << listeners.size() << endl;
}

void SurfData::removeListener(Surface* surface)
{
  //// We had some problems if we were removing things from the list
  //// while it was being iterated over
  ////if (!cleanupStarted) {
  //  // make sure its OK to erase the object, then do so
  //  list<Surface*>::iterator itr =
  //    find(listeners.begin(),listeners.end(),surface);
  //  //cout << "remove: " << surface << endl;
  //  if(itr!=listeners.end()) {
  //    //cout << "remove: " << *itr << endl;
  //    //listeners.erase(itr);
  //    *itr = 0;
  //    cout << "Listener removed: " << listeners.size()
  //         << "this: " << this << endl;
  //  }
  //  for (list<Surface*>::iterator itr2 = listeners.begin();
  //        itr2 != listeners.end();
  //        ++itr2) {
  //    cout << "Still listening: " << *itr << endl;
  //  }
  ////} else {
  ////  cout << "Cleanup has already started" << endl;
  ////}
  cout << "SurfData: " << this
       << " #listbef: " << listeners.size()
       << " removing: " << surface;
  listeners.remove(surface);
  cout << " #listaft: " << listeners.size() << endl;
}

/// For use with copy constructor and assignment operator-- creates a list of
/// pointers to the points in the data set which is used to check for 
/// duplication when other points are added in the future
void SurfData::buildOrderedPoints()
{
  orderedPoints.clear();
  for (unsigned i = 0; i < points.size(); i++) {
    orderedPoints.insert(points[i]);
  }
}

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
      xMatrix[pt+xval*mapping.size()] = points[mapping[pt]]->X()[xval];
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
    throw surfpack::file_open_failure(filename);
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
    throw surfpack::file_open_failure(filename);
  } else if (binary) {
    readBinary(infile);
  } else {
    readText(infile);
  }
  // Object may have already been created
  valid.xMatrix = false;
  valid.yVector = false;
  infile.close();
}

/// Write a set of SurfPoints to an output stream
void SurfData::writeBinary(ostream& os) const 
{
  //try {
    unsigned s = mapping.size();
    os.write((char*)&s,sizeof(s));
    os.write((char*)&xsize,sizeof(xsize));
    os.write((char*)&fsize,sizeof(fsize));
    for (unsigned i = 0; i < mapping.size(); i++) {
      points[mapping[i]]->writeBinary(os);
    }
  //} catch(...) {
  //  cerr << "Exception caught and rethrown in SurfData::writeBinary" << endl;
  //  throw;
  //}
}

/// Write a set of SurfPoints to an output stream
void SurfData::writeText(ostream& os) const
{
  //try {
    // Write an extra space after last token to prevent
    // ifstream from setting the bad_bit when data is later read in.
    //os << mapping.size() << " " << endl
    //   << xsize << " " << endl 
    //   << fsize << " " << endl;
    os << mapping.size() << endl
       << xsize << endl 
       << fsize << endl;
    for (unsigned i = 0; i < mapping.size(); i++) {
      points[mapping[i]]->writeText(os);
    }
  //} catch(...) {
  //  cerr << "Exception caught and rethrown in SurfData::writeText" << endl;
  //  throw;
  //}
}

/// Read a set of SurfPoints from an input stream
void SurfData::readBinary(istream& is) 
{
  unsigned numPointsRead = 0;
  unsigned size;
  try {
    cleanup();
    is.read((char*)&size,sizeof(size));
    is.read((char*)&xsize,sizeof(xsize));
    is.read((char*)&fsize,sizeof(fsize));
    points.clear();
    for (numPointsRead = 0; numPointsRead < size; numPointsRead++) {
      // True for second argument signals a binary read
      surfpack::checkForEOF(is);
      //points.push_back(SurfPoint(xsize,fsize,is,true));  
      this->addPoint(SurfPoint(xsize,fsize,is,true));  
    }
    defaultMapping();
  } catch(surfpack::io_exception&) {
    cerr << "Expected: " << size << " points.  "
         << "Read: " << numPointsRead << " points." << endl;
    throw;
  } //catch(...) {
  //  cerr << "Exception caught and rethrown in SurfData::readBinary" << endl;
  //  throw;
  //}
}

/// Read a set of SurfPoints from an input stream
void SurfData::readText(istream& is) 
{
  unsigned numPointsRead = 0;
  unsigned size;
  try {
    cleanup();
    string sline;
    getline(is,sline);
    istringstream streamline(sline);
    streamline >> size;
    getline(is,sline);
    streamline.str(sline);
    streamline.clear();
    streamline >> xsize;
    getline(is,sline);
    streamline.str(sline);
    streamline.clear();
    streamline >> fsize;
    points.clear();
    for (numPointsRead = 0; numPointsRead < size; numPointsRead++) {
      surfpack::checkForEOF(is);
      // False for last argument signals a text read
      //points.push_back(SurfPoint(xsize,fsize,is,false));  
      this->addPoint(SurfPoint(xsize,fsize,is,false));  
    }
    defaultMapping();
  } catch(surfpack::io_exception&) {
    cerr << "Expected: " << size << " points.  "
         << "Read: " << numPointsRead << " points." << endl;
    throw;
  } //catch(...) {
  //  cerr << "Exception caught and rethrown in SurfData::readText" << endl;
  //  throw;
  //}
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
  if (surfpack::hasExtension(filename,".sd")) {
    return true;
  } else if (surfpack::hasExtension(filename,".txt")) {
    return false;
  } else {
    throw surfpack::io_exception(
      "Unrecognized filename extension.  Use .sd or .txt"
    );
  }
}

void SurfData::notifyListeners(int msg)
{
  //cout << "Size: " << listeners.size() << endl;
  if (listeners.size() != 0) {
    cout << "SurfData: " << this
         << " thinks it has " << listeners.size()
         << " surfaces to notify: " << endl;
  }
  list<Surface*>::iterator itr = listeners.begin();
  while (itr != listeners.end()) {
    cout << "\tnotifying: " << *itr << " of " << msg << endl;
    if (*itr) {
      (*itr)->notify(msg);
    }
    ++itr;
  }
}

// ____________________________________________________________________________
// Testing 
// ____________________________________________________________________________

// Throw an exception if there are any mismatches in the number of
// dimensions or number of response values among points in the data set  
void SurfData::sanityCheck() const
{
  if (!points.empty()) {
    unsigned dimensionality = points[0]->xSize();
    unsigned numResponses = points[0]->fSize();
    for (unsigned i = 1; i < points.size(); i++) {
      if (points[i]->xSize() != dimensionality ||
          points[i]->fSize() != numResponses ) {
        ostringstream errormsg;
        errormsg << "Error in SurfData::sanityCheck." << endl
		 << "Point 0 has " << dimensionality << " dimensions "
                 << "and " << numResponses << " response values, " << endl
                 << "but point " << i << " has " << points[i]->xSize()
 		 << " dimensions and " << points[i]->fSize() << "response "
 		 << " values.";
	throw bad_surf_data(errormsg.str());
      } // if mismatch
    } // for each point
  } // if !empty
}

/// Check that the index falls within acceptable boundaries (i.e., is
/// less than mapping.size()
void SurfData::checkRangeNumPoints(const string& header, unsigned index) const
{
  if (index >= mapping.size()) {
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

void SurfData::checkRangeNumResponses(const string& header, unsigned index) const
{
  if (index >= fsize) {
    ostringstream errormsg;
    errormsg << header << endl;
    if (fsize == 0) {
      errormsg << "Index " << index << " specified, but there are zero response"
	       << "values."
               << endl;
    } else {
      errormsg << "Requested: " 
	     << index 
	     << "; actual max index: "
	     << fsize - 1
	     << endl;
    }
    throw range_error(errormsg.str());
  }
}

#ifdef __TESTING_MODE__
  int SurfData::constructCount = 0;
  int SurfData::copyCount = 0;
  int SurfData::destructCount = 0;
#endif

