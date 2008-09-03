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
#include "ModelFactory.h"
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

SurfpackModelFactory* ModelFactory::createModelFactory(ParamMap& args)
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
