// These lines were automatically prepended
#ifndef SURFPACK_CONFIG_H
#define SURFPACK_CONFIG_H
#include "config.h"
#endif
// End of prepended lines

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
#include <set>

#include "SurfPoint.h"
#include "SurfData.h"
#include "surfpack.h"
#include "AxesBounds.h"



using namespace std;


AxesBounds::AxesBounds(vector<AxesBounds::Axis> axes_) 
  : axes(axes_), ndims(axes.size()), point(ndims)
{

}

AxesBounds::AxesBounds(string filename) 
//void readInput(istream& is = cin)
{
    ifstream infile(filename.c_str(), ios::in);
    if (!infile) {
      throw surfpack::file_open_failure(filename);
    } 
    npts = 1;
    string sline;
    getline(infile, sline);
    istringstream istream(sline);
    istream >> ndims;
    axes.resize(ndims);
    point.resize(ndims);
    for (int i = 0; i < ndims; i++) {
        getline(infile, sline);
	istringstream streamline(sline);
        char c;
        streamline >> c;
	if (c == 'v') {
	    streamline >> axes[i].min >> axes[i].max >> axes[i].pts;
            if (axes[i].pts > 1) {
              axes[i].interval = (axes[i].max - axes[i].min) / (axes[i].pts - 1);
            } else {
              axes[i].interval = 0.0;
            }
	    npts *= axes[i].pts;
        } else if (c == 'f') {
            streamline >> axes[i].min;
	    axes[i].max = axes[i].min;
	    axes[i].pts = 1;
	    axes[i].interval = 0.0;
	} else {
    	    cerr << "Expected 'f' or 'v'" << endl;
	}
    }
    infile.close();
}    
     
AxesBounds::~AxesBounds()
{

}

void AxesBounds::initialize()
{
    point.resize(ndims);
    for (int i = 0; i < ndims; i++) {
	point[i] = 0;
    }
}

void AxesBounds::nextPoint()
{
    int curDim = ndims-1;
    while (axes[curDim].pts == 1 || point[curDim] == (axes[curDim].pts - 1)) {
	curDim--;
    }
    if (curDim < ndims) {
        point[curDim]++;
	curDim++;
	while(curDim < ndims) {
	    point[curDim] = 0;
	    curDim++;
	}
    }
}

//SurfData* populate()
//{
//    initialize();
//    surfptx.resize(ndims);
//    for (int i = 0; i < npts; i++) {
//        for (int j = 0; j < ndims; j++) {
//	    surfptx[j] = (axes[j].max - axes[j].min)*((double)rand()/(double)INT_MAX)+ axes[j].min;
//    	    //cout << setw(10) << surfptx[j];
//	}
//	//cout << endl;
//	//double fx = surface ? surface->evaluate(surfptx) : rastrigin(surfptx)*.5;
//        SurfPoint sp(surfptx);
//	surfData->addPoint(sp);
//	nextPoint();
//    }
//    //cout << "Number of points: " << npts << endl;
//}

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
  return  new SurfData(sps);
}

SurfData* AxesBounds::sampleMonteCarlo(unsigned numPts, const vector<string>& testFunctions) 
{
  surfptx.resize(ndims);
  vector<SurfPoint> sps;
  for (unsigned i = 0; i < numPts; i++) {
      for (unsigned j = 0; j < ndims; j++) {
	surfptx[j] = (axes[j].max - axes[j].min)*((double)rand()/(double)INT_MAX)+ axes[j].min;
      }
      SurfPoint sp(surfptx);
      for (unsigned k = 0; k < testFunctions.size(); k++) {
        sp.addResponse(surfpack::testFunction(testFunctions[k], sp.X()));
      }
      sps.push_back(sp);
      //surfData->addPoint(sp);
      nextPoint();
  }
  return new SurfData(sps);
}


//void randomSample(vector< string >& args) 
//{
//  SurfData sd = randomPoints(args[1]);
//  if (args.size() == 4) {
//    vector<double> newResponseValues(sd.size());
//    //unsigned newindex = sd->addResponse();
//    for (unsigned i = 0; i < sd.size(); i++) {
//      newResponseValues[i] = surfpack::testFunction(args[3],sd[i].X());
//      //sd->Point(i).F(newindex,response);
//    }
//    sd.addResponse(newResponseValues);
//  }
//  sd.write(args[2]);
//}

//void gridPoints(vector< string >& args) 
//{
//  SurfData sd = SurfDataGenerator::pointSpecToSurfData(args[1]);
//  if (args.size() == 4) {
//    vector<double> newResponseValues(sd.size());
//    //unsigned newindex = sd->addResponse();
//    for (unsigned i = 0; i < sd.size(); i++) {
//      newResponseValues[i] = surfpack::testFunction(args[3],sd[i].X());
//      //double response = testFunction(args[3],sd->Point(i).X());
//      //sd->Point(i).F(newindex,response);
//    }
//    sd.addResponse(newResponseValues);
//  }
//  sd.write(args[2]);
//}

