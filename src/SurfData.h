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

class SurfPoint;
//class Surface;

class SurfData
{
// ____________________________________________________________________________
// Creation, Destruction, Initialization 
// ____________________________________________________________________________

public:

  /// Vector of points will be copied
  SurfData(const std::vector<SurfPoint>& points);

  /// Read a set of SurfPoints from a file
  SurfData(const std::string filename);

  /// Makes a deep copy of the object 
  SurfData(const SurfData& sd); 
  
  /// STL data members' resources automatically deallocated 
  ~SurfData();
  
private:

  // Default constructor explicitly disallowed. 
  SurfData();

// ____________________________________________________________________________
// Overloaded operators 
// ____________________________________________________________________________
public:

  /// makes deep copy 
  SurfData& operator=(const SurfData& sd);

  /// makes deep comparison
  bool operator==(const SurfData& sd);

  /// makes deep comparison
  bool operator!=(const SurfData& sd);

// ____________________________________________________________________________
// Queries 
// ____________________________________________________________________________

  /// Return the number of SurfPoints in the data set 
  unsigned size() const;
  
  /// Return the dimensionality of the SurfPoints 
  unsigned xSize() const;

  /// Return the number of response values in the data set
  unsigned fSize() const;

  /// Return a point from the data set
  SurfPoint& Point(unsigned index);

  /// Return a reference to the SurfPoints vector 
  std::vector<SurfPoint>& Points();

// ____________________________________________________________________________
// Commands 
// ____________________________________________________________________________

  /// Add a point to the data set. The parameter point will be copied.
  void addPoint(const SurfPoint& sp);

  /// Add a new response variable to each point. 
  /// Return the index of the new variable.
  unsigned addResponse(); 
  
  /// Inform this object that a Surface wants to be notified when this object
  /// changes
  //void addListener(Surface *);
 
  /// remove the Surface from the list of surfaces that are notified when the
  /// data changes
  //void removeListener(Surface *);

   
// ____________________________________________________________________________
// I/O
// ____________________________________________________________________________

  /// Write a set of SurfPoints to a file
  void write(const std::string filename) const;

  /// Write a set of SurfPoints to an output stream
  std::ostream& write(std::ostream& os = std::cout, bool binary = false) const;

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

  /// List of pointers to listening/observing Surface objects 
  /// which need to be notified when this object changes
  //std::list<Surface*> listeners;

// ____________________________________________________________________________
// Helper methods 
// ____________________________________________________________________________

  /// Notify listening surfaces whenever a SurfPoint
  /// is added or removed.
  //void notifyListeners(); 

// ____________________________________________________________________________
// Testing 
// ____________________________________________________________________________

#ifdef __TESTING_MODE__ 
  friend class SurfDataUnitTest;
  friend class SurfaceUnitTest;

public:
  static int constructCount;
  static int copyCount;
  static int destructCount;
#endif
};

/// so a SurfData object can be printed
std::ostream& operator<<(std::ostream& os, const SurfData& data);

#endif
