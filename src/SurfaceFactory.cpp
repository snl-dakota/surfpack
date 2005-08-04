#include "surfpack_config.h"
#include "surfpack.h"
#include "SurfaceFactory.h"
#include "ANNSurface.h"
#include "KrigingCPPSurface.h"
#include "KrigingSurface.h"
#include "MarsSurface.h"
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

/// The parameter type specifies the name of a Surface class.  Instantiate
/// the appropriate class using parameter surfData, and return the new
/// surface.  The client is responsible to call delete on the Surface* that
/// is returned.
Surface* SurfaceFactory::createSurface(const string& type, SurfData* sd)
{
  if (type == "Polynomial") {
    return new PolynomialSurface(sd);
  } else if (type == "Kriging") {
    return new KrigingSurface(sd);
  } else if (type == "KrigingCPP") {
    return new KrigingCPPSurface(sd);
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

