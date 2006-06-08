/*  _______________________________________________________________________

    DAKOTA: Design Analysis Kit for Optimization and Terascale Applications
    Copyright (c) 2001, Sandia National Laboratories.

    Surfpack: A Software Library of Multidimensional Surface Fitting Methods

    Surfpack is distributed under the DAKOTA GNU General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */
#include "surfpack_config.h"
#include "surfpack.h"
#include "SurfaceFactory.h"
#include "ANNSurface.h"
#include "KrigingCPPSurface.h"
#include "KrigingSurface.h"
#include "MarsSurface.h"
#include "MarsValidator.h"
#include "PolynomialSurface.h"
#include "RBFNetSurface.h"

class SurfData;
using namespace std;

/// Open up the file specified by parameter filename.  The first item in the 
/// file should be the name of the surface type.  Once that information is 
/// known, build create and return the appropriate surface.  The client is
/// responsible to call delete on the Surface* that is returned.
Surface* SurfaceFactory::createSurface(const string& filename)
{
  const string name = surfpack::surfaceName(filename);
  if (name == "polynomial") {
    return new PolynomialSurface(filename); 
  } else if (name == "kriging") {
    return new KrigingSurface(filename);
  } else if (name == "mars") {
    return new MarsSurface(filename);
  } else if (name == "marsc") {
    return new MarsCppSurface(filename);
  } else if (name == "rbf") {
    return new RBFNetSurface(filename);
  } else if (name == "ann") {
    return new ANNSurface(filename);
  } else {
    cerr << "Unknown surface type: " << name << endl;
    return 0;
  }
}

/// The parameter type specifies the name of a Surface class.  Instantiate
/// the appropriate class using parameter surfData, and return the new
/// surface.  The client is responsible to call delete on the Surface* that
/// is returned.
Surface* SurfaceFactory::createSurface(const string& type, SurfData* sd)
{
  if (type == "polynomial") {
    return new PolynomialSurface(sd);
  } else if (type == "kriging") {
    return new KrigingSurface(sd);
  } else if (type == "kriging_cpp") {
    return new KrigingCPPSurface(sd);
  } else if (type == "mars") {
    return new MarsSurface(sd);
  } else if (type == "marsc") {
    return new MarsCppSurface(sd);
  } else if (type == "rbf") {
    return new RBFNetSurface(sd);
  } else if (type == "ann") {
    return new ANNSurface(sd);
  } else {
    ostringstream os;
    os << "Unknown surface type: " << type; 
    throw os.str(); 
  }
}

