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
SurfData::SurfData(const SurfData& sd) : xsize(sd.xsize), fsize(sd.fsize),
  points(sd.points), excludedPoints(sd.excludedPoints), 
  mapping(sd.mapping), valid(sd.valid)
{
#ifdef __TESTING_MODE__
  constructCount++;
  copyCount++;
#endif
  if (valid.xMatrix) {
    unsigned numElements = mapping.size() * xsize;
    xMatrix = new double[numElements];
    memcpy(xMatrix,sd.xMatrix, numElements*sizeof(double));
  } else {
    xMatrix = 0;
  }
 
  if (valid.yVector) {
    unsigned numElements = mapping.size();
    yVector = new double[numElements];
    memcpy(yVector,sd.yVector, numElements*sizeof(double));
  } else {
    yVector = 0;
  }
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
  delete xMatrix;
  delete yVector;
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
  newSD->setDefaultIndex(defaultIndex);
  return newSD;
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

    if (valid.xMatrix) {
      unsigned numElements = mapping.size() * xsize;
      xMatrix = new double[numElements];
      memcpy(xMatrix,sd.xMatrix, numElements*sizeof(double));
    } else {
      xMatrix = 0;
    }
 
    if (valid.yVector) {
      unsigned numElements = mapping.size();
      yVector = new double[numElements];
      memcpy(yVector,sd.yVector, numElements*sizeof(double));
    } else {
      yVector = 0;
    }

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
const SurfPoint& SurfData::operator[](unsigned index) const
{
  checkRange(index);
  return points[mapping[index]];
}

/// Return a point from the data set
//SurfPoint& SurfData::Point(unsigned index) 
SurfPoint& SurfData::operator[](unsigned index) 
{
  checkRange(index);
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
std::set<unsigned> SurfData::getExcludedPoints() const 
{
  return excludedPoints;
}

/// Get the response value of the (index)th point that corresponds to this
/// surface
double SurfData::getResponse(unsigned index) const
{
  checkRange(index);
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

// ____________________________________________________________________________
// Commands 
// ____________________________________________________________________________

/// Specify which response value getResponse will return
void SurfData::setDefaultIndex(unsigned index)
{
  if (index >= fsize) {
    cerr << "Cannot set defaultIndex; " << index << " exceeds # of response" 
         << endl;
  } else {
    // If the defaultIndex is being changed, it invalidates the yVector
    if (defaultIndex != index) {
      valid.yVector = false;
    }
    defaultIndex = index;
  }
}
  
/// Set the response value of the (index)th point that corresponds to this
/// surface
void SurfData::setResponse(unsigned index, double value)
{
  checkRange(index);
  if (points[mapping[index]].F(defaultIndex) != value) {
    points[mapping[index]].F(defaultIndex, value);
    valid.yVector = false;
  }
}
  
/// Add a point to the data set. The parameter point will be copied.
void SurfData::addPoint(const SurfPoint& sp) 
{
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
  if (points.empty()) {
    cerr << "Cannot add another response because there are no points "
         << "in the data set." << endl;
    return 0;
  } else if (newValues.size() != points.size()) {
    cerr << "Cannot add another response because the number of response"
         << " values does not match the size of the physical data set." 
         << " All physical data points must be included in order to add"
         << " a response variable" 
         << endl;
  } else {
    newindex = points[0].addResponse(newValues[0]);
    fsize = newindex + 1;
    for (unsigned i = 1; i < points.size(); i++) {
      if (points[i].addResponse(newValues[i]) != newindex) {
        cerr << "Size mismatch among data points in this SurfData object."
             << endl
             << "This will likely cause a fatal error." 
             << endl;
      }
    }
  }
  return newindex;
}

/// Specify which points should be skipped
void SurfData::setExcludedPoints(std::set<unsigned> excludedPoints)
{
  if (excludedPoints.empty()) {
    defaultMapping();
  } else {
    // Remove any values from the excluded points set that are out of
    // range
    for (set<unsigned>::iterator itr = excludedPoints.begin();
         itr != excludedPoints.end();
         ++itr) {
      if (*itr >= points.size()) {
        cerr << "Excluded point index exceeds set size" << endl;
        excludedPoints.erase(itr);
      } 
    }
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
    if (mappingIndex != mapping.size()) {
      cerr << "Miscalculated mapping" << endl;
    }
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

/// Make sure an index falls within acceptable boundaries
void SurfData::checkRange(unsigned index) const
{
  // What if there are no points?
  if (points.empty()) {
    cerr << "Cannot return points[" << index << "].  There are no points." << endl;
  } else if (index >= mapping.size()) {
    cerr << "Out of range in SurfData::checkRange" << endl;
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
void SurfData::write(const std::string filename) const
{
  bool binary = (filename.find(".txt") != filename.size() - 4);
  ofstream outfile(filename.c_str(), (binary ? ios::out|ios::binary : ios::out));
  if (!outfile) {
    cerr << "File named \"" 
         << filename
	 << "\" could not be opened." 
	 << endl;
    return;
  } else if (binary) {
    writeBinary(outfile);
  } else {
    writeText(outfile);
  }
  outfile.close();
}

/// Read a set of SurfPoints from a file.  Opens file and calls other version.
void SurfData::read(const string filename)
{
  bool binary = (filename.find(".txt") != filename.size() - 4);
  ifstream infile(filename.c_str(), (binary ? ios::in|ios::binary : ios::in));
  if (!infile) {
    cerr << "File named \"" 
         << filename
	 << "\" could not be opened." 
	 << endl;
    return;
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
  unsigned s = mapping.size();
  os.write((char*)&s,sizeof(s));
  os.write((char*)&xsize,sizeof(xsize));
  os.write((char*)&fsize,sizeof(fsize));
  for (unsigned i = 0; i < mapping.size(); i++) {
    points[mapping[i]].writeBinary(os);
  }
}

/// Write a set of SurfPoints to an output stream
void SurfData::writeText(ostream& os) const
{
  // Using trailing words to circumvent istringstream bug
  //os << points.size() << " data points" << endl
  //   << xsize << " input variables" << endl 
  //   << fsize << " response variables" << endl;
  os << mapping.size() << " " << endl
     << xsize << " " << endl 
     << fsize << " " << endl;
  for (unsigned i = 0; i < mapping.size(); i++) {
    points[mapping[i]].writeText(os);
  }
}

/// Read a set of SurfPoints from an input stream
void SurfData::readBinary(istream& is) 
{
  unsigned size;  
  is.read((char*)&size,sizeof(size));
  is.read((char*)&xsize,sizeof(xsize));
  is.read((char*)&fsize,sizeof(fsize));
  points.clear();
  for (unsigned i = 0; i < size; i++) {
    // True for second argument signals a binary read
    points.push_back(SurfPoint(xsize,fsize,is,true));  
  }
  defaultMapping();
}

/// Read a set of SurfPoints from an input stream
void SurfData::readText(istream& is) 
{
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
    // False for second argument signals a text read
    points.push_back(SurfPoint(xsize,fsize,is,false));  
  }
  defaultMapping();
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

//void SurfData::notifyListeners()
//{
//  list<Surface*>::iterator itr;
//  for(itr=listeners.begin();itr!=listeners.end();++itr)
//    (*itr)->notify();
//}

// ____________________________________________________________________________
// Testing 
// ____________________________________________________________________________

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

