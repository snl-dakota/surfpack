#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include "SurfData.h"
#include "Surface.h"
#include "PolynomialSurface.h"
#include "MarsSurface.h"
#include "KrigingSurface.h"
#include "RBFNetSurface.h"
#include "ANNSurface.h"
#include "surfpack.h"

//______________________________________________________________________________
// Functions for use by Surface methods
//______________________________________________________________________________

using namespace std;

const string surfaceName(const string filename)
{
  bool binary = (filename.find(".txt") != filename.size() - 4);
  ifstream infile(filename.c_str(), (binary ? ios::in|ios::binary : ios::in));
  if (!infile) {
    cerr << "Could not open " << filename << "." << endl;
    return "none";
  } else if (binary) {
    unsigned nameSize;
    infile.read(reinterpret_cast<char*>(&nameSize),sizeof(nameSize));
    char* surfaceType = new char[nameSize+1];
    infile.read(surfaceType,nameSize);
    surfaceType[nameSize] = '\0';
    infile.close();
    return string(surfaceType);
  } else {
    string nameInFile;
    getline(infile,nameInFile);
    infile.close();
    return nameInFile;
  }
    
} 

Surface* createSurface(const string filename)
{
  const string name = surfaceName(filename);
  if (name == "Polynomial") {
    return new PolynomialSurface(filename); 
  } else if (name == "Kriging") {
    return new KrigingSurface(filename);
  } else if (name == "Mars") {
    return new MarsSurface(filename);
  } else if (name == "RBFNet") {
    return new RBFNetSurface(filename);
  } else if (name == "ANN") {
    return new ANNSurface(filename);
  } else {
    cerr << "Unknown surface type: " << name << endl;
    return 0;
  }
}

Surface* createSurface(const std::string type, 
  AbstractSurfDataIterator* dataItr)
{
  if (type == "Kriging") {
    return new KrigingSurface(dataItr);
  } else if (type == "Mars") {
    return new MarsSurface(dataItr);
  } else if (type == "RBFNet") {
    return new RBFNetSurface(dataItr);
  } else if (type == "ANN") {
    return new ANNSurface(dataItr);
  } else {
    cerr << "Unknown surface type: " << type << endl;
    return 0;
  }
}

Surface* createSurface(const std::string type,
  AbstractSurfDataIterator* dataItr, unsigned order)
{
  if (type == "Polynomial") {
    return new PolynomialSurface(dataItr, order); 
  } else {
    cerr << "Unknown surface type: " << type << endl;
    return 0;
  }
}

Surface* createSurface(const std::string type, SurfData& sd, 
  unsigned responseIndex)
{
  if (type == "Kriging") {
    return new KrigingSurface(sd, responseIndex);
  } else if (type == "Mars") {
    return new MarsSurface(sd, responseIndex);
  } else if (type == "RBFNet") {
    return new RBFNetSurface(sd, responseIndex);
  } else if (type == "ANN") {
    return new ANNSurface(sd, responseIndex);
  } else {
    cerr << "Unknown surface type: " << type << endl;
    return 0;
  }
}

Surface* createSurface(const std::string type, SurfData& sd, 
  unsigned responseIndex, unsigned order)
{

  if (type == "Polynomial") {
    return new PolynomialSurface(sd, order, responseIndex); 
  } else {
    cerr << "Unknown surface type: " << type << endl;
    return 0;
  }
}

double euclideanDistance(const vector<double>& pt1, const vector<double>& pt2)
{
  double distance = 0.0;
  if (pt1.size() != pt2.size()) {
    cerr << "Cannot compute euclidean distance.  Vectors have different sizes."
         << endl;
  } else {
    for (unsigned i = 0; i < pt1.size(); i++) {
      distance += (pt1[i] - pt2[i]) * (pt1[i] - pt2[i]);
    }
    distance = sqrt(distance);
  }
  return distance;

}

void vectorDifference(vector<double>& diff, const vector<double>& pt1,
  const vector<double>& pt2)
{
  if (pt1.size() != pt2.size() || pt1.size() != diff.size()) {
    cerr << "Cannot compute vector difference: size mismatch." << endl;
    return;
  }
  for (unsigned i = 0; i < pt1.size(); i++) {
    diff[i] = pt1[i] - pt2[i];
  }
}


void printVector(const std::string header, vector<double>& vec)
{
  cout << header << " size: " << vec.size() << endl;
  for (unsigned i = 0; i < vec.size(); i++) {
    cout << i << " " << vec[i] << endl;
  }
}
