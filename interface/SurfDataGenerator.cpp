/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#ifdef HAVE_CONFIG_H
#include "surfpack_config.h"
#endif
#include "surfpack.h"
#include "SurfData.h"
#include "SurfDataGenerator.h"

using std::cerr;
using std::endl;
using std::ifstream;
using std::ios;
using std::istream;
using std::istringstream;
using std::ostream;
using std::setw;
using std::string;
using std::vector;

void SurfDataGenerator::readPointSpec(vector<Axis>& axes, istream& is)
{
    // read the number of dimensions
    string line;
    getline(is, line);
    istringstream istream(line);
    int ndims;
    istream >> ndims;
    axes.resize(ndims);

    // read specification for each dimension
    for (int i = 0; i < axes.size(); i++) {
        getline(is, line);
	istringstream streamline(line);
        char c;
        streamline >> c;
	if (c == 'v') {
            // TO DO: What happens if last value isn't there?
            // TO DO: What happens if min >= max ?
	    streamline >> axes[i].min >> axes[i].max >> axes[i].pts;
            axes[i].interval = (axes[i].max - axes[i].min) / (axes[i].pts - 1);
        } else if (c == 'f') {
            streamline >> axes[i].min;
	    axes[i].max = axes[i].min;
	    axes[i].pts = 1;
	    axes[i].interval = 0.0;
	} else {
    	    cerr << "Expected 'f' or 'v'" << endl;
	}
    }
}    

SurfData SurfDataGenerator::pointSpecToSurfData(string filename) {
  vector<double> surfptx;
  ifstream infile(filename.c_str(),ios::in);
  if (!infile) {
    throw surfpack::file_open_failure(filename);
  }
  vector<Axis> axes;
  readPointSpec(axes,infile);
  infile.close();
  //SurfDataGenerator::printPointSpec(axes,cout);
  PointCounter pc(axes);
  pc.initialize();
  surfptx.resize(axes.size());
  vector<SurfPoint> sps;
  while (!pc.atLastPoint()) {
      for (int j = 0; j < axes.size(); j++) {
          surfptx[j] = axes[j].min + axes[j].interval*pc[j];
  	    //cout << setw(10) << surfptx[j];
      }
      //cout << endl;
      //double fx = surface ? surface->evaluate(surfptx) : rastrigin(surfptx)*.5;
      SurfPoint sp(surfptx);
      sps.push_back(sp);
      //surfData->addPoint(sp);
      //pc.print(cout);
      pc.nextPoint();
  }
  return SurfData(sps);
}

void SurfDataGenerator::pointSpecToSurfDataFile(string pointSpecFilename,
  string outputFilename, string functionName)
{
  vector<string> container(1);
  container[0] = functionName;
  pointSpecToSurfDataFile(pointSpecFilename, outputFilename, container);
}

void SurfDataGenerator::pointSpecToSurfDataFile(string pointSpecFilename,
  string outputFilename, vector<string> functionNames)
{
  SurfData sd = pointSpecToSurfData(pointSpecFilename);
  vector<double> responses(sd.size());
  for (unsigned i = 0; i < functionNames.size(); i++) {
    for (unsigned j = 0; j < responses.size(); j++) {
      responses[j] = surfpack::testFunction(functionNames[i], sd[j].X());
    }
    sd.addResponse(responses);
  }
  sd.write(outputFilename);
}

void SurfDataGenerator::pointSpecToSurfDataFile(string pointSpecFilename,
  string outputFilename)
{
  SurfData sd = pointSpecToSurfData(pointSpecFilename);
  sd.write(outputFilename); 
}

void SurfDataGenerator::printPointSpec(vector<Axis>& axes, ostream& os)
{
  os << axes.size() << endl;
  for (unsigned i = 0; i < axes.size(); i++) {
    Axis& a = axes[i];
    if (a.min == a.max) {
      os << "f " << a.min << endl;
    } else {
      os << "v " 
	 << setw(8) << a.min 
	 << setw(8) << a.max 
	 << setw(8) << a.pts
	 << setw(8) << a.interval
	 << endl;
    }
  } 

}
     
SurfDataGenerator::PointCounter::PointCounter(vector<Axis>& axes_) : axes(axes_)
{
  initialize();
}

SurfDataGenerator::PointCounter::~PointCounter()
{

}

void SurfDataGenerator::PointCounter::initialize()
{
  point.resize(axes.size());
  for (int i = 0; i < axes.size(); i++) {
    point[i] = 0;
  }
  lastPoint = false;
}

void SurfDataGenerator::PointCounter::nextPoint()
{
  //unsigned nDimsMaxedOut = 0;
  if (!lastPoint) {
    int curDim = axes.size()-1;
    while (point[curDim] == (axes[curDim].pts - 1) && curDim >= 0) {
      curDim--;
    }
    // When the last dimension needs to "roll over" to zero, you're at the last point
    if (curDim < 0) {
      lastPoint = true;
      return;
    }
    if (point[curDim] < axes[curDim].pts - 1) {
      point[curDim]++;
      //if (point[curDim] == axes[curDim].pts - 1) {
      //  nDimsMaxedOut++;
      //}
    }
    curDim++;
    while(curDim < axes.size()) {
      point[curDim] = 0;
      //nDimsMaxedOut = 0;
      curDim++;
    }
    //if (nDimsMaxedOut = axes.size()) {
    //  lastPoint = true;
    //}
  }
}

void SurfDataGenerator::PointCounter::print(ostream& os) const
{
  for (unsigned i = 0; i < axes.size(); i++) {
    os << setw(3) << point[i];
  }
  os << endl;
} 

bool SurfDataGenerator::PointCounter::atLastPoint() const
{
  return lastPoint;
}

const double& SurfDataGenerator::PointCounter::operator[](unsigned i) const
{
  return point[i];
}

//SurfData randomPoints(string filename) {
//  ifstream infile(filename.c_str(),ios::in);
//  if (!infile) {
//          cerr << "File " << filename << " not found." << endl;
//          exit(1);
//  }
//  readInput(infile);
//  infile.close();
//  initialize();
//  surfptx.resize(axes.size());
//  vector<SurfPoint> sps;
//  for (int i = 0; i < npts; i++) {
//      for (int j = 0; j < axes.size(); j++) {
//	surfptx[j] = (axes[j].max - axes[j].min)*((double)rand()/(double)INT_MAX)+ axes[j].min;
//  	    //cout << setw(10) << surfptx[j];
//      }
//      //cout << endl;
//      //double fx = surface ? surface->evaluate(surfptx) : rastrigin(surfptx)*.5;
//      SurfPoint sp(surfptx);
//      sps.push_back(sp);
//      //surfData->addPoint(sp);
//      nextPoint();
//  }
//  return SurfData(sps);
//}
//
//
//double testFunction(const string name, const vector<double>& pt)
//{
//  if (name == "rosenbrock") {
//    return surfpack::rosenbrock(pt);
//  } else if (name == "sphere") {
//    return surfpack::sphere(pt);
//  } else {
//    return surfpack::rastrigin(pt);
//  }
//}

