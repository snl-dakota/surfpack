// These lines were automatically prepended
#ifndef SURFPACK_CONFIG_H
#define SURFPACK_CONFIG_H
#include "config.h"
#endif
// End of prepended lines

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <set>
#include <cmath>
#include "SurfData.h"
#include "Surface.h"
#include "PolynomialSurface.h"
#include "MarsSurface.h"
#include "KrigingSurface.h"
#include "RBFNetSurface.h"
#include "ANNSurface.h"
#include "surfpack.h"
#include "SurfaceFactory.h"

//______________________________________________________________________________
// Functions for use by Surface methods
//______________________________________________________________________________

using namespace std;

Surface* SurfaceFactory::createSurface(const string& filename)
{
  const string name = surfpack::surfaceName(filename);
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

Surface* SurfaceFactory::createSurface(const string& type, SurfData* sd)
{
  if (type == "Polynomial") {
    return new PolynomialSurface(sd);
  } else if (type == "Kriging") {
    return new KrigingSurface(sd);
  } else if (type == "Mars") {
    return new MarsSurface(sd);
  } else if (type == "RBFNet") {
    return new RBFNetSurface(sd);
  } else if (type == "ANN") {
    return new ANNSurface(sd);
  } else {
    ostringstream os;
    os << "Unknown surface type: " << type; 
    throw os.str(); 
  }
}

//Surface* SurfaceFactory::createSurface(const std::string*, SurfData* sd, 
//    const SurfpackParser::ArgList& arglist)
//{ 
//  //Surface* surf = createSurface(type, sd);
//  //surf->config(arglist);
//  //return surf;
//  return 0;
//}

//Surface* SurfaceFactory::createSurface(const string& type, SurfData& sd, unsigned order)
//{
//  if (type == "Polynomial") {
//    return new PolynomialSurface(&sd, order); 
//  } else {
//    cerr << "Unknown surface type: " << type << endl;
//    return 0;
//  }
//}

