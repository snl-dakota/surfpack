/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */
#ifndef AXES_BOUNDS_H
#define AXES_BOUNDS_H

#include "surfpack_config.h"
#include "surfpack_system_headers.h"

class SurfData;
class SurfPoint;

/// Concrete class used in conjunction with Surfpack commands GridPoints and 
/// MonteCarloSample.  Minimum and maximum values are specified along each
/// dimension.  A SurfData object can be created by specifying some number of
/// simple Monte Carlo samples from hypercube defined by the boundaries.
/// Alternatively, a client may specify a number of points along each dimension
/// (in addition to the maximum and minimum) and then create a SurfData object 
/// from that hypergrid of points.
class AxesBounds
{
public:
  /// Values for one dimension For random samples, only bounds are used.  For
  /// hyper-gridding, the pts and interval fields are also used.  
  struct Axis {
      /// Minimum value along this dimension.
      double min;

      /// Maximum value along this dimension.
      double max;

      /// The number of raster points along the dimension (used for 
      /// hyper-gridding).
      int pts;

      /// The spacing needed between max and min to fit pts points along this 
      /// dimension.
      double interval;
  };

  enum ParamType {file = 1, data = 2};

  /// Object is created from an existing list of Axis objects
  AxesBounds( std::vector<Axis> axes_);

  /// Object is created from a text file
  AxesBounds( std::string fileOrData, ParamType pt = file);

  /// No special behavior
  ~AxesBounds();

  /// Return a hypergrid data set as a SurfData object.  The client is 
  /// responsible to deallocate the memory.
  SurfData* sampleGrid(const std::vector<std::string>& testFunctions);

  /// Return a data set with numPts SurfPoints.  Parameter testFunctions
  /// must contain the names of zero or more functions at which all the data
  /// points should be evaluated.  The test functions should reside in the
  /// surfpack namespace.
  SurfData* sampleMonteCarlo(unsigned numPts, 
    const std::vector<std::string>& testFunctions);

  /// Reset the counter used to iterate through dimensions for grid data.
  /// The next point to be added to the data set would be (axes[0].min, 
  /// axes[1].min, ..., axes[ndims - 1].min).
  void initialize();

  /// Advance the counter used to iterate through dimensions for grid data.
  /// For example, if the client has requested a 10 x 10 grid, and the point
  /// data member currently holds (3,9), it will be advanced (not unlike an
  /// odometer) to (4,0).  This will signify that the next point to be added
  /// should be (axes[0].min+4*axes[0].interval, axes[1].min+0*axes[1].interval.
  void nextPoint();

protected:
  void parseBounds(std::istream& is);

protected:
  /// Counter used to iterate through the dimensions for grid data.  The 
  /// sampleGrid method uses it in conjunction with the min, max, and interval
  /// values along each axis to create a SurfData object.  For example, if
  /// a value of (3,2) for the point data member would signify that the next
  /// point to be added to the grid should be (axes[0].min+3*axes[0].interval,
  /// axes[1].min+2*axes[1].interval
  std::vector<int> point;

  /// Used during the iteration that creates a SurfData set.  It holds the
  /// location of the current SurfPoint that is to be added.
  std::vector<double> surfptx;
  
  /// The set of <minimum, maximum, #intervals> specifications for each
  /// dimension
  std::vector<Axis> axes;

  /// Number of dimensions in the data set
  unsigned ndims;

  /// Number of points to be in the SurfData object returned by sampleGrid
  unsigned npts;
};
  
#endif
