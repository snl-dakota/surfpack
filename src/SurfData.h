/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#ifndef __SURF_DATA_H__
#define __SURF_DATA_H__

#ifdef HAVE_CONFIG_H
#include "surfpack_config.h"
#endif
#include "surfpack_system_headers.h"

#include "SurfPoint.h"
#include "SurfpackParser.h"

class Surface;
class SurfScaler;

/// Contains a set of SurfPoint objects.  May be associated with zero or more
/// Surface objects, which it notifies when its data changes or when it goes 
/// out of existence.  Contains support for exclusion of some of the SurfPoint
/// objects that are physically present, so that clients may operate on some
/// subset of the data if they so choose (e.g., in a cross-validation 
/// algorithm).  Contains methods for I/O support.  Does not allow duplicate
/// points.
/// \todo Allow the points to be weighted differently, as would be needed
/// in weighted regression.
class SurfData
{
public:
/// Nested exception class. A bad_surf_data exception is thrown whenever a 
/// client attempts to:
/// 1) Add a SurfPoint to the SurfData object that has a different number of
///    dimensions and/or a different number of response values than another
///    SurfPoint already in the set;
/// 2) invoke addResponse(...) when there are not yet any SurfPoints;
/// 3) invoke addResponse(...) with the wrong number of values;
/// 4) invoke addResponse(...) on an object where the logical size of the 
///    data set does not match the phyiscal size (i.e., some of the SurfPoints
///    have been marked for exclusion);
/// 5) write a SurfData object that contains no data to a stream;
/// 6) create a Surface with a SurfData object that does not have enough points. 
class bad_surf_data : public std::runtime_error
{
public:
  bad_surf_data(const std::string& msg = "") : std::runtime_error(msg) {}
};

// ____________________________________________________________________________
// Creation, Destruction, Initialization 
// ____________________________________________________________________________

public:
  /// Vector of points will be copied and checked for duplicates
  SurfData(const std::vector<SurfPoint>& points_);

  /// Read a set of SurfPoints from a file
  SurfData(const std::string filename);

  /// Read a set of SurfPoints from a std::istream
  SurfData(std::istream& is, bool binary = false);

  /// Read a set of SurfPoints from a std::istream.  The stream does not
  /// contain the normal header information (#points, #vars, #responses).
  /// The #vars and #responses are explicitly specified in the constructor;
  /// The stream reader processes data until eof, assuming one point per line.
  SurfData(const std::string filename, unsigned n_vars, unsigned n_responses, 
    unsigned n_cols_to_skip);

  /// Make a deep copy of the object 
  SurfData(const SurfData& other); 
  
  /// STL data members' resources automatically deallocated 
  ~SurfData();

  /// Copy only the points which have not been marked for exclusion
  SurfData copyActive();
  
  /// First SurfPoint added will determine the dimensions of the data set 
  SurfData();

private:
  /// Data member initialization that is common to all constructors
  void init();
  
  /// Call delete on the SurfPoint* in the data set.
  void cleanup();

// ____________________________________________________________________________
// Overloaded operators 
// ____________________________________________________________________________
public:
  /// Makes deep copy 
  SurfData& operator=(const SurfData& other);

  /// Makes deep comparison
  bool operator==(const SurfData& other) const;

  /// Makes deep comparison
  bool operator!=(const SurfData& other) const;

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

  /// Return the number of response functions in the data set
  unsigned fSize() const;

  /// Returns true if the data has been scaled
  bool isScaled() const;

  /// Return the set of excluded points (the indices)
  const std::set<unsigned>& getExcludedPoints() const ; 

  /// Get the response value of the (index)th point
  double getResponse(unsigned index) const;

  /// Return defaultIndex
  unsigned getDefaultIndex() const;

  /// Retrieve the label for one of the predictor variables
  const std::string& getXLabel(unsigned index) const;

  /// Retrieve the label for one of the predictor variables
  const std::string& getFLabel(unsigned index) const;

  /// Retrieve the index and variable type (predictor/response) for a given
  /// name.  Return false if not found
  bool varIndex(const std::string& name, unsigned& index, bool& isResponse) const;

  

// ____________________________________________________________________________
// Commands 
// ____________________________________________________________________________

  /// Specify which response value getResponse will return. When a Surface 
  /// object that is associated with the SurfData object operates on the data,
  /// it sets this value so that the response value lookup function will return
  /// the value for the response variable that that particular Surface object
  /// is interested in.  
  void setDefaultIndex(unsigned index); 
  
  /// Set the response value of the (index)th point that corresponds to this
  /// surface
  void setResponse(unsigned index, double value);

  /// Calculates parameters so that the data can be viewed as scaled
  void setScaler(SurfScaler* scaler_in);
  
  /// Add a point to the data set. The parameter point will be copied.
  void addPoint(const SurfPoint& sp);

  /// Add a new response variable to each point. 
  /// Return the index of the new variable.
  unsigned addResponse(const std::vector<double>& newValues, 
    std::string label = ""); 
  
  /// Specify which points should be skipped.  This can be used when only a 
  /// subset of the SurfPoints should be used for some computation.
  void setExcludedPoints(const std::set<unsigned>& excluded_points);

  /// Inform this object that a Surface wants to be notified when this object
  /// changes
  void addListener(Surface*);
 
  /// Remove the Surface from the list of surfaces that are notified when the
  /// data changes
  void removeListener(Surface*);

  /// For use with copy constructor and assignment operator-- creates a list of
  /// pointers to the points in the data set which is used to check for 
  /// duplication when other points are added in the future
  void buildOrderedPoints();

  /// Iterate through the points, setting the current scaler
  void enableScaling();

  /// Iterate through the points, clearing the scaler
  void disableScaling();

  /// Set the labels for the predictor variables
  void setXLabels(const std::vector<std::string>& labels);

  /// Set the labels for the response variables
  void setFLabels(const std::vector<std::string>& labels);

  /// Set the label for a single response variable
  void setFLabel(unsigned index, const std::string& label);

private:
  /// Maps all indices to themselves in the mapping data member
  void defaultMapping();

  /// Set x vars labels to 'x0' 'x1', etc.; resp. vars to 'f0' 'f1', etc.
  void defaultLabels();

public:
   
// ____________________________________________________________________________
// I/O
// ____________________________________________________________________________

  /// Write a set of SurfPoints to a file
  void write(const std::string& filename) const;

  /// Read a set of SurfPoints from a file
  void read(const std::string& filename);
  
  /// Write the surface in binary format
  void writeBinary(std::ostream& os) const;

  /// Write the surface in text format
  void writeText(std::ostream& os, bool write_header = true,
    bool write_labels = true) const ; 

  /// Read the surface in binary format
  void readBinary(std::istream& is); 

  /// Read the surface in text format
  void readText(std::istream& is, bool read_header = true, 
    unsigned skip_columns = 0); 

private:

// ____________________________________________________________________________
// Data members 
// ____________________________________________________________________________

  /// Dimensionality of the space from wich the SurfPoints are drawn
  unsigned xsize;

  /// Number of response variables in the data set 
  unsigned fsize;

  /// Controls how the data is scaled during processing
  SurfScaler* scaler;

  /// The set of points in this data set
  std::vector<SurfPoint*> points; 

  /// The indices of points that are to be excluded in computation. This can
  /// be used in a cross-validation scheme to systematically ignore parts of
  /// data set at different times.  
  std::set<unsigned> excludedPoints;

  /// For mapping the indices in points to the indices returned by operator[].
  /// Normally, mapping[i] is equal to i, but the set of excludedPoints is not
  /// empty, this will not be the case.
  std::vector<unsigned> mapping;

  /// The index of the response variable that will be returned by F
  mutable unsigned defaultIndex;

  /// Labels for the predictor variables
  std::vector< std::string > xLabels;
 
  /// Labels for the response variables
  std::vector< std::string > fLabels;

public:
  typedef std::set<SurfPoint*,SurfPoint::SurfPointPtrLessThan> SurfPointSet;

private:
  /// Stores the same set of SurfPoint* that points does, but because it is a 
  /// set, membership tests can be done in O(log n) time.  When combined with
  /// the SurfPointPtrLessThan functor object, it allows a SurfData object to
  /// check all SurfPoints in the data set against all others for duplication.
  /// This can be done in O(n log n) time instead of O(n^2).
  SurfPointSet orderedPoints;

  /// List of pointers to listening/observing Surface objects that need to be 
  /// notified when this object changes
  std::list<Surface*> listeners;

// ____________________________________________________________________________
// Constants 
// ____________________________________________________________________________
public:
  /// Used to send a message through Surface::notify(...) that this object is
  /// going out of existence
  static const int GOING_OUT_OF_EXISTENCE;

  /// Used to send a message through Surface::notify(...) that one or more
  /// SurfPoints have been added or modified
  static const int DATA_MODIFIED;
  
// ____________________________________________________________________________
// Helper methods 
// ____________________________________________________________________________

  /// Notify listening surfaces whenever something of interest happens to this
  /// data set
  void notifyListeners(int msg); 

  /// Returns true if file has .bspd extension, false if it has .spd extension. 
  /// Otherwise, an exception is thrown.
  bool hasBinaryFileExtension(const std::string& filename) const;

  /// If the line contains single-quoted string, parse them out as labels
  /// and return true; otherwise, return false
  bool readLabelsIfPresent(std::string single_line);

  /// Read the #points, #vars, #responses
  unsigned SurfData::readHeaderInfo(std::istream& is);

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

#ifdef __TESTING_MODE__ 
  friend class SurfDataTest;
  friend class SurfaceTest;
#endif
};

/// Print the SurfData to a stream 
std::ostream& operator<<(std::ostream& os, const SurfData& data);

#endif
