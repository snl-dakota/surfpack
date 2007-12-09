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
#include "SurfaceFactory.h"
#include "ANNSurface.h"
#include "KrigingSurface.h"
#include "MarsSurface.h"
#include "PolynomialSurface.h"
#include "RBFNetSurface.h"

#include "SurfpackModel.h"
#include "LinearRegressionModel.h"
#include "RadialBasisFunctionModel.h"
#include "DirectANNModel.h"
#include "KrigingModel.h"
#include "MovingLeastSquaresModel.h"
#include "MarsModel.h"

class SurfData;
using std::cerr;
using std::endl;
using std::ostringstream;
using std::string;

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
  } else if (type == "mars") {
    return new MarsSurface(sd);
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

SurfpackModelFactory* SurfaceFactory::createModelFactory(ParamMap& args)
{
  string type = args["type"];
  SurfpackModelFactory* smf; 
  if (type == "") throw string("Model must declare type");
  else if (type == "polynomial") {
    smf = new LinearRegressionModelFactory(args);
  } else if (type == "mls") {
    smf = new MovingLeastSquaresModelFactory(args);
  } else if (type == "rbf") {
    smf = new RadialBasisFunctionModelFactory(args);
  } else if (type == "kriging") {
    smf = new KrigingModelFactory(args);
  } else if (type == "ann") {
    smf = new DirectANNModelFactory(args);
  } else if (type == "mars") {
    smf = new MarsModelFactory(args);
  } else {
    throw string("Model type requested not recognized");
  }
  return smf;
}
