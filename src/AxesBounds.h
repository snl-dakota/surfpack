// These lines were automatically prepended
#ifndef SURFPACK_CONFIG_H
#define SURFPACK_CONFIG_H
#include "config.h"
#endif
// End of prepended lines

#ifndef POINT_DEFINITION_H
#define POINT_DEFINITION_H

class SurfData;
class SurfPoint;

#include <string>
#include <vector>

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
  /// Values for one dimension.  The pts field is used to specify the number
  /// of points along the dimension.  The interval field can be used to compute
  /// the increment that must be used to get the right number of points between
  /// the minimum and maximum values.
  struct Axis {
      double min;
      double max;
      int pts;
      double interval;
  };

  /// Object is created from an existing list of Axis objects
  AxesBounds( std::vector<Axis> axes_);

  /// Object is created from a text file
  AxesBounds( std::string filename);

  /// No special behavior
  ~AxesBounds();
  /// Return a hypergrid data set.  Client is responsible to deallocate the
  /// memory.
  SurfData* sampleGrid(const std::vector<std::string>& testFunctions);

  /// Return a data set with numPts SurfPoints.  Parameter testFunctions
  /// must contain the names of zero or more functions at which all the data
  /// points should be evaluated.  The test functions should reside in the
  /// surfpack namespace.
  SurfData* sampleMonteCarlo(unsigned numPts, 
    const std::vector<std::string>& testFunctions);

  /// Reset the counter used to iterate through dimensions for grid data
  void initialize();

  /// Advance the counter used to iterate through dimensions for grid data
  void nextPoint();

protected:
  /// Counter used to iterate through the dimensions for grid data.  The 
  /// sampleGrid method uses it in conjunction with the min, max, and interval
  /// values along each axis to create a SurfData object.
  std::vector<int> point;
  std::vector<double> surfptx;
  std::vector<Axis> axes;
  unsigned ndims;
  unsigned npts;
};
  
#endif
