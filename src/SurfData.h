// These lines were automatically prepended
#ifndef SURFPACK_CONFIG_H
#define SURFPACK_CONFIG_H
#include "config.h"
#endif
// End of prepended lines

// ----------------------------------------------------------
// Project: SURFPACK++
//
// File:        SurfData.h
// Author:      Eric Cyr
// Modified:	Mark Richards
//
// Description: 
// + SurfData class - this is the container for all the data 
//   from which an empirical model is created. 
// + Left shift (<<) operator for SurfData. 
// ----------------------------------------------------------

#ifndef __SURF_DATA_H__
#define __SURF_DATA_H__

#include <stdexcept>

class SurfPoint;
//class Surface;

class SurfData
{
private:
class bad_surf_data : public std::runtime_error
{
public:
  bad_surf_data(const std::string& msg = "") : std::runtime_error(msg) {}
};

struct SurfDataStateConsistency 
{
  bool xMatrix;
  bool yVector;
  SurfDataStateConsistency() : xMatrix(false), yVector(false) {}
};

// ____________________________________________________________________________
// Creation, Destruction, Initialization 
// ____________________________________________________________________________

public:
  /// Vector of points will be copied
  SurfData(const std::vector<SurfPoint>& points);

  /// Read a set of SurfPoints from a file
  SurfData(const std::string filename);

  /// Read a set of SurfPoints from an istream
  SurfData(std::istream& is, bool binary = false);

  /// Makes a deep copy of the object 
  SurfData(const SurfData& other); 
  
  /// STL data members' resources automatically deallocated 
  ~SurfData();

  /// Copy only the "active" points
  SurfData copyActive();
  
protected:
  /// Default constructor explicitly disallowed by clients, but subclassing 
  /// permitted.
  SurfData();

private:
  /// Initialize data members
  void init();
  
  /// Copy xMatrix and yVector from another SurfData object
  void copyBlockData(const SurfData& other);

  /// Deallocate any memory allocated for xMatrix and/or yVector
  void cleanup();

// ____________________________________________________________________________
// Overloaded operators 
// ____________________________________________________________________________
public:
  /// makes deep copy 
  SurfData& operator=(const SurfData& sd);

  /// makes deep comparison
  bool operator==(const SurfData& sd) const;

  /// makes deep comparison
  bool operator!=(const SurfData& sd) const;

  /// Return a const reference to SurfPoint at given index
  const SurfPoint& operator[](unsigned index) const;

// ____________________________________________________________________________
// Queries 
// ____________________________________________________________________________

  /// Return the number of SurfPoints in the data set 
  unsigned size() const;

  /// True if there are no points
  bool empty() const;
  
  /// Return the dimensionality of the SurfPoints 
  unsigned xSize() const;

  /// Return the number of response values in the data set
  unsigned fSize() const;

  /// Return a point from the data set
  //SurfPoint& Point(unsigned index);

  /// Return a reference to the SurfPoints vector 
  //std::vector<SurfPoint>& Points();

  /// Return the set of excluded points (the indices)
  const std::set<unsigned>& getExcludedPoints() const ; 

  /// Get the response value of the (index)th point that corresponds to this
  /// surface
  double getResponse(unsigned index) const;

  /// Get default index
  unsigned getDefaultIndex() const;

  /// Return point domains as a matrix in a contiguous block.  Be careful.
  /// The data should not be changed.
  const double* getXMatrix() const;

  /// Return response values for the default response in a contiguous block.  
  /// Be careful. The data should not be changed.
  const double* getYVector() const;

  /// Returns true if the filename extension is .sd.
  static bool hasBinaryExtension(const std::string& filename);
  
  /// Returns true if the filename extension is .txt.
  static bool hasTextExtension(const std::string& filename);
  
  

// ____________________________________________________________________________
// Commands 
// ____________________________________________________________________________

  /// Specify which response value getResponse will return
  void setDefaultIndex(unsigned index); 
  
  /// Set the response value of the (index)th point that corresponds to this
  /// surface
  void setResponse(unsigned index, double value);
  
  /// Add a point to the data set. The parameter point will be copied.
  void addPoint(const SurfPoint& sp);

  /// Add a new response variable to each point. 
  /// Return the index of the new variable.
  unsigned addResponse(const std::vector<double>& newValues); 
  
  /// Specify which points should be skipped
  void setExcludedPoints(const std::set<unsigned>& excludedPoints);

  /// Inform this object that a Surface wants to be notified when this object
  /// changes
  //void addListener(Surface *);
 
  /// remove the Surface from the list of surfaces that are notified when the
  /// data changes
  //void removeListener(Surface *);
private:
  /// Maps all indices to themselves
  void defaultMapping();

  /// Creates a matrix of the domains for all of the points in a contiguous
  /// block of memory, for use in matrix operations
  void validateXMatrix() const;
 
  /// Creates a vector of response values for the default response value in
  /// a contiguous blcok of memory
  void validateYVector() const;
public:
   
// ____________________________________________________________________________
// I/O
// ____________________________________________________________________________

  /// Write a set of SurfPoints to a file
  void write(const std::string& filename) const;

  /// Write a set of SurfPoints to an output stream
  void read(const std::string& filename);
  
  /// Write the surface in binary format
  void writeBinary(std::ostream& os) const;

  /// Write the surface in text format
  void writeText(std::ostream& os) const ; 

  /// Read the surface in binary format
  void readBinary(std::istream& is); 

  /// Read the surface in text format
  void readText(std::istream& is); 

private:

// ____________________________________________________________________________
// Data members 
// ____________________________________________________________________________

  /// Dimensionality of the space from wich the SurfPoints are drawn
  unsigned xsize;

  /// Number of response variables in the data set 
  unsigned fsize;

  /// The set of points in this data set
  std::vector<SurfPoint> points; 

  /// The indices of points in points that are skipped
  std::set<unsigned> excludedPoints;

  /// For mapping the indices in points to the indices returned by operator[]
  std::vector<unsigned> mapping;

  /// Pointer to the domain of the data points, represented as a contiguous 
  /// block of memory
  mutable double* xMatrix;

  /// Pointer to a set of response values for the data poitns, stored in a
  /// contiguous block of memory
  mutable double* yVector;

  /// The index of the response variable that will be returned by F
  mutable unsigned defaultIndex;

  /// Keeps track of which data members are valid
  mutable SurfDataStateConsistency valid;



  /// List of pointers to listening/observing Surface objects 
  /// which need to be notified when this object changes
  //std::list<Surface*> listeners;

// ____________________________________________________________________________
// Helper methods 
// ____________________________________________________________________________

  /// Notify listening surfaces whenever a SurfPoint
  /// is added or removed.
  //void notifyListeners(); 

 /// Returns true if file is opened in binary mode, false if it is opened
 /// in text mode.  If file cannot be opened, an exception is thrown.
 bool testFileExtension(const std::string& filename) const;

// ____________________________________________________________________________
// Testing 
// ____________________________________________________________________________

protected:
  // Throw an exception if there are any mismatches in the number of
  // dimensions or number of response values among points in the data set  
  void sanityCheck() const;

  /// Make sure an index falls within acceptable boundaries
  void checkRangeNumPoints(const std::string& header, unsigned index) const;

  /// Make sure an index falls within acceptable boundaries
  void checkRangeNumResponses(const std::string& header, unsigned index) const;

public:

#ifdef __TESTING_MODE__ 
  friend class SurfDataTest;
  friend class SurfaceTest;

  static int constructCount;
  static int copyCount;
  static int destructCount;
#endif
};

/// so a SurfData object can be printed
std::ostream& operator<<(std::ostream& os, const SurfData& data);

#endif
