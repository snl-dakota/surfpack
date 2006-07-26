/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#include "surfpack_config.h"

#include "surfpack.h"
#include "SurfPoint.h"
#include "SurfData.h"
#include "AxesBounds.h"

using namespace std;
using surfpack::dbg;
const int axdbg = 0;

/// Object is created from an existing list of Axis objects
AxesBounds::AxesBounds(vector<AxesBounds::Axis> axes_in) 
  : axes(axes_in)
{

}

/// Object is created from a text file
AxesBounds::AxesBounds(string file_or_data, ParamType pt) 
{
    // Consider a file with the following data:
    // 0 1 |
    // 2 |
    // -1 1

    // Each line of the file corresponds to one variable
    // If there are two values on the line, they represent the minimum
    // and maximum value for that variable.  If there is only one value
    // that means that the variable is fixed (fixed variable = oxymoron?).
    // White space is actually ignored.  The ranges for all dimensions
    // could be one the same line.  The vertical bars between the values
    // for different dimensions are required though.

    // Make sure file is opened successfully.
    // Note: infile is only used inside case file, but it cannot
    // be declared inside a switch statement
    assert(pt == file || pt == data);
    ifstream infile(file_or_data.c_str(), ios::in);
    switch(pt) {
      case file: // the string is a file name
        if (!infile) {
          throw surfpack::file_open_failure(file_or_data);
        } 
        parseBounds(infile);
        break;
      case data: // the string actually holds the boundary data
        istringstream is(file_or_data);
        parseBounds(is);
        break;
    }
    infile.close();
}    

/// Read (min,max) pairs or a fixed value for each dimension.
/// Data for different dimensions should be separated by '|' 
void AxesBounds::parseBounds(std::istream& is)
{
  axes.push_back(Axis());
  string token;
  while (!is.eof()) {
    // read in min for current dimension
    is >> axes.back().min;
    // read in next token -- it should be either a '|' or the max for this dim
    is >> token;
    dbg(axdbg) << "Token read; " << token << " eof?: " << is.eof() << "\n";
    if (is.eof()) break;
    if (token == "|") {
      axes.back().max = axes.back().min;
      axes.push_back(Axis());
      continue; // There is no 'max' for this dim; skip to next one
    } else {
      axes.back().max = atof(token.c_str());
      axes.back().minIsMax = false;
      // now the next token should be a '|' or eof
      is >> token;
    dbg(axdbg) << "Token read; " << token << " eof?: " << is.eof() << "\n";
      if (is.eof()) break;
      if (token != "|") {
        cerr << "Expected |" << endl;
 	exit(1);
      }
      axes.push_back(Axis());
    }
  }

  if (axdbg) { // debug output
    cout << "Axes values parsed" << endl;
    for (unsigned i = 0; i < axes.size(); i++) {
      cout << axes[i].min;
      if (!axes[i].minIsMax) {
        cout << " " << axes[i].max;
      }
      cout << endl;
    } 
    cout << "dims: " << axes.size() << endl;
  } // debug output
      
}

/// Based on the ranges for each variable contained in axes and the number
/// of grid points requested per dimension (grid_points), computes and returns t/// the interval between two values on each dimension.
vector<double> AxesBounds::computeIntervals(vector<Axis>& axes, 
  const vector<unsigned>& grid_points)
{
  assert(axes.size() == grid_points.size());
  vector<double> intervals(axes.size());
  for (unsigned i = 0; i < grid_points.size(); i++) {
    // Treating one as special case avoids div. by zero below
    if (grid_points[i] == 1) {
      intervals[i] = 0.0;
    } else {
      dbg(axdbg) << "i " << i << " min/max: " << axes[i].min << " " << axes[i].max << " gp: " << grid_points[i] << " int: ";
      intervals[i] = (axes[i].max - axes[i].min)/(grid_points[i] - 1);
      dbg(axdbg) << intervals[i] << "\n";
    }
  }
  return intervals;
}     

/// No special behavior
AxesBounds::~AxesBounds()
{

}


/// Advance the counter used to iterate through dimensions for grid data.
/// For example, if the client has requested a 10 x 10 grid, and the odometer 
/// data member currently holds (3,9), it will be advanced (not unlike an
/// odometer) to (4,0).  This will signify that the next point to be added
/// should be (axes[0].min+4*intervals[0], axes[1].min+0*intervals[0].
void AxesBounds::nextPoint(vector<unsigned>& point_odometer, 
  const vector<unsigned>& grid_points)
{
    // Scan across the "odometer" reading to find the first value that does
    // not need to roll over.  For example, if the user requested a 5 x 5 x 5 x
    // 5 grid and point holds (3,2,4,4), then the next value should be 
    // (3,3,0,0), so the rightmost two values need to roll over, and cur_dim,
    // should point to the 2.
    int cur_dim = axes.size()-1;
    while (grid_points[cur_dim] == 1 || 
	   point_odometer[cur_dim] == (grid_points[cur_dim] - 1)) {
	cur_dim--;
    }

    // If the odometer isn't maxed out, increase the digit at cur_dim by one,
    // and then zero out everything to the right.
    if (cur_dim < axes.size()) {
        point_odometer[cur_dim]++;
	cur_dim++;
	while(cur_dim < axes.size()) {
	    point_odometer[cur_dim] = 0;
	    cur_dim++;
	}
    }
}

/// Return a hypergrid data set as a SurfData object.  The client is 
/// responsible to deallocate the memory.
SurfData* AxesBounds::sampleGrid(const vector<unsigned>& grid_points, const vector<string>& test_functions) 
{
  vector<unsigned> point_odometer(grid_points.size(),0);
  vector<double> surfptx(axes.size());
  vector<SurfPoint> sps;
  vector<double> intervals = computeIntervals(axes,grid_points);
  unsigned npts = accumulate(grid_points.begin(),grid_points.end(),
	1, multiplies<unsigned>());
  for (int i = 0; i < npts; i++) {
      for (int j = 0; j < axes.size(); j++) {
          surfptx[j] = axes[j].min + intervals[j]*point_odometer[j];
      }
      SurfPoint sp(surfptx);
      for (unsigned k = 0; k < test_functions.size(); k++) {
        sp.addResponse(surfpack::testFunction(test_functions[k], sp.X()));
      }
      sps.push_back(sp);
      nextPoint(point_odometer, grid_points);
  }
  SurfData* sd = new SurfData(sps);
  if (!test_functions.empty()) sd->setFLabels(test_functions);
  return sd;
}

/// Return a data set with size SurfPoints.  Parameter test_functions
/// must contain the names of zero or more functions at which all the data
/// points should be evaluated.  The test functions should reside in the
/// surfpack namespace.
SurfData* AxesBounds::sampleMonteCarlo(unsigned size, 
  const vector<string>& test_functions) 
{
  vector<double> surfptx(axes.size());
  vector<SurfPoint> sps;
  for (unsigned i = 0; i < size; i++) {
      for (unsigned j = 0; j < axes.size(); j++) {
	surfptx[j] = (axes[j].max - axes[j].min) * 
          ((double)rand()/(double)INT_MAX)+ axes[j].min;
      }
      SurfPoint sp(surfptx);
      for (unsigned k = 0; k < test_functions.size(); k++) {
        sp.addResponse(surfpack::testFunction(test_functions[k], sp.X()));
      }
      sps.push_back(sp);
  }
  return new SurfData(sps);
}
