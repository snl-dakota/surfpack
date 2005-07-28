#include "config.h"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "surfpack.h"
#include "SurfPoint.h"
#include "SurfData.h"
#include "AxesBounds.h"

using namespace std;

/// Object is created from an existing list of Axis objects
AxesBounds::AxesBounds(vector<AxesBounds::Axis> axes_) 
  : axes(axes_), ndims(axes.size()), point(ndims)
{

}

/// Object is created from a text file
AxesBounds::AxesBounds(string fileOrData, ParamType pt) 
{
    // Consider a file with the following data:

    // 3
    // v 0 10 10
    // v 0 10 10
    // f 5

    // The first line indicates that the data set will have three dimensions.
    // The next line says that the grid should go from 0 to 10 along the first
    // dimension, with 10 data points.  Because the end points will always be
    // included, the interval between points will be a little over 1.  (If it
    // were exactly one, there would be points at 0, 1,..., 10, which would be
    // 11 points total.  The second dimension also has ten points.  Since these
    // are the only "variable" dimensions in the data set, the SurfData object
    // will have 10 x 10 = 100 points.  The last line says that the third 
    // dimension is fixed at the value 5.  So all 100 data points will have
    // x[2] = 3.
    
    // Make sure file is opened successfully.
    assert(pt == file || pt == data);
    ifstream infile(fileOrData.c_str(), ios::in);
    switch(pt) {
      case file:
        if (!infile) {
          throw surfpack::file_open_failure(fileOrData);
        } 
        parseBounds(infile);
        break;
      case data:
        replace(fileOrData.begin(),fileOrData.end(),'|','\n');
        istringstream is(fileOrData);
        parseBounds(is);
        break;
    }
    infile.close();
}    

void AxesBounds::parseBounds(std::istream& is)
{

    // Read the number of dimensions.
    npts = 1;
    string sline;
    getline(is, sline);
    istringstream istream(sline);
    istream >> ndims;
    axes.resize(ndims);
    point.resize(ndims);
    // Read the values for each dimension
    for (int i = 0; i < ndims; i++) {
        getline(is, sline);
	istringstream streamline(sline);
        char c;
        streamline >> c;
        // Is this a variable or fixed dimension?
	if (c == 'v') {
            // If it's variable, read min, max, and pts.
            // Then compute the value for interval.
	    streamline >> axes[i].min >> axes[i].max >> axes[i].pts;
            if (axes[i].pts > 1) {
              axes[i].interval = (axes[i].max - axes[i].min)/(axes[i].pts - 1);
            } else {
              axes[i].interval = 0.0;
            }
	    npts *= axes[i].pts;
        } else if (c == 'f') {
            // If it's fixed, read only min and max.  Use dummy values of 1
            // and 0.0 for the pts and interval fields.
            streamline >> axes[i].min;
	    axes[i].max = axes[i].min;
	    axes[i].pts = 1;
	    axes[i].interval = 0.0;
	} else {
            // Each line after the first MUST begin with 'f' or 'v'.
            ostringstream msg;
            msg << "Expected 'f' or 'v' on line " 
                << i+1
                << "in AxesBounds::parseBounds" 
                << endl;
 	    throw surfpack::io_exception(msg.str());
	}
    }
}
     
/// No special behavior
AxesBounds::~AxesBounds()
{

}

/// Reset the counter used to iterate through dimensions for grid data.
/// The next point to be added to the data set would be (axes[0].min, 
/// axes[1].min, ..., axes[ndims - 1].min).
void AxesBounds::initialize()
{
    point.resize(ndims);
    for (int i = 0; i < ndims; i++) {
	point[i] = 0;
    }
}

/// Advance the counter used to iterate through dimensions for grid data.
/// For example, if the client has requested a 10 x 10 grid, and the point
/// data member currently holds (3,9), it will be advanced (not unlike an
/// odometer) to (4,0).  This will signify that the next point to be added
/// should be (axes[0].min+4*axes[0].interval, axes[1].min+0*axes[1].interval.
void AxesBounds::nextPoint()
{
    // Scan across the "odometer" reading to find the first value that does
    // not need to roll over.  For example, if the user requested a 5 x 5 x 5 x
    // 5 grid and point holds (3,2,4,4), then the next value should be 
    // (3,3,0,0), so the rightmost two values need to roll over, and curDim,
    // should point to the 2.
    int curDim = ndims-1;
    while (axes[curDim].pts == 1 || point[curDim] == (axes[curDim].pts - 1)) {
	curDim--;
    }

    // If the odometer isn't maxed out, increase the digit at curDim by one,
    // and then zero out everything to the right.
    if (curDim < ndims) {
        point[curDim]++;
	curDim++;
	while(curDim < ndims) {
	    point[curDim] = 0;
	    curDim++;
	}
    }
}

/// Return a hypergrid data set as a SurfData object.  The client is 
/// responsible to deallocate the memory.
SurfData* AxesBounds::sampleGrid(const vector<string>& testFunctions) 
{
  initialize();
  surfptx.resize(ndims);
  vector<SurfPoint> sps;
  for (int i = 0; i < npts; i++) {
      for (int j = 0; j < ndims; j++) {
          surfptx[j] = axes[j].min + axes[j].interval*point[j];
      }
      SurfPoint sp(surfptx);
      for (unsigned k = 0; k < testFunctions.size(); k++) {
        sp.addResponse(surfpack::testFunction(testFunctions[k], sp.X()));
      }
      sps.push_back(sp);
      nextPoint();
  }
  return new SurfData(sps);
}

/// Return a data set with numPts SurfPoints.  Parameter testFunctions
/// must contain the names of zero or more functions at which all the data
/// points should be evaluated.  The test functions should reside in the
/// surfpack namespace.
SurfData* AxesBounds::sampleMonteCarlo(unsigned numPts, 
  const vector<string>& testFunctions) 
{
  surfptx.resize(ndims);
  vector<SurfPoint> sps;
  for (unsigned i = 0; i < numPts; i++) {
      for (unsigned j = 0; j < ndims; j++) {
	surfptx[j] = (axes[j].max - axes[j].min) * 
          ((double)rand()/(double)INT_MAX)+ axes[j].min;
      }
      SurfPoint sp(surfptx);
      for (unsigned k = 0; k < testFunctions.size(); k++) {
        sp.addResponse(surfpack::testFunction(testFunctions[k], sp.X()));
      }
      sps.push_back(sp);
      nextPoint();
  }
  return new SurfData(sps);
}
