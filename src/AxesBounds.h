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

/// Concrete class used in conjunction with Surfpack command CreateSample. 
/// Minimum and maximum values are specified along each dimension in the 
/// parameter space.  A SurfData object can be created by requesting
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
      Axis() : min(0.0), max(0.0), minIsMax(true) {}

      /// Minimum value along this dimension.
      double min;

      /// Maximum value along this dimension.
      double max;

      /// No variation along this dimension
      bool minIsMax;
  };

  enum ParamType {file = 1, data = 2};

  /// Object is created from an existing list of Axis objects
  AxesBounds( std::vector<Axis> axes_in);

  /// Object is created from a text file
  AxesBounds( std::string file_or_data, ParamType pt = file);

  /// No special behavior
  ~AxesBounds();

  /// Return a hypergrid data set as a SurfData object.  The client is 
  /// responsible to deallocate the memory.
  SurfData* sampleGrid(const std::vector<unsigned>& grid_points, 
    const std::vector<std::string>& test_functions);

  /// Return a data set with size SurfPoints.  Parameter test_functions
  /// must contain the names of zero or more functions at which all the data
  /// points should be evaluated.  The test functions should reside in the
  /// surfpack namespace.  The client must deallocate the memory.
  SurfData* sampleMonteCarlo(unsigned size, 
    const std::vector<std::string>& test_functions);

  /// Advance the counter used to iterate through dimensions for grid data.
  /// For example, if the client has requested a 10 x 10 grid, and the point
  /// data member currently holds (3,9), it will be advanced (not unlike an
  /// odometer) to (4,0).  This will signify that the next point to be added
  /// should be (axes[0].min+4*axes[0].interval, axes[1].min+0*axes[1].interval.
  void nextPoint(std::vector<unsigned>& point_odometer,
    const std::vector<unsigned>& grid_points);

  std::vector<double> computeIntervals(std::vector<Axis>& axes, 
    const std::vector<unsigned>& grid_points);

protected:
  /// Parse out the boundary information from an input stream created from a 
  /// string or file
  void parseBounds(std::istream& is);

protected:
  /// The set of <minimum, maximum> specifications for each dimension
  std::vector<Axis> axes;
};
#endif
