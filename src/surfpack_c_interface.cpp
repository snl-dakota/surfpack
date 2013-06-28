/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#include "surfpack.h"
#include "surfpack_c_interface.h"
#include "SurfpackModel.h"

// Must include the headers of any models that might be loaded, or
// they won't be registered with the serializer
#include "DirectANNModel.h"
#include "KrigingModel.h"
#include "LinearRegressionModel.h"
#include "MarsModel.h"
#include "MovingLeastSquaresModel.h"
#include "RadialBasisFunctionModel.h"

/* Implementation of simplified C interface to Surfpack libraries */


/// one global instance of a SurfpackModel
static SurfpackModel* surfpackCModel = NULL;

// Model load example copied from SurfpackInterpreter.cpp; 
extern "C" 
int surfpack_load_model(const char * const model_filename)
{
  bool success = true;
  try {

    bool binary = surfpack::isBinaryModelFilename(model_filename);

    std::ifstream model_ifstream(model_filename);
    if (!model_ifstream.good())
      throw "Failure opening model file for load."; 

    if (binary) {
      std::cout << "Loading model binary file '" << model_filename << "'." 
		<< std::endl;
      boost::archive::binary_iarchive input_archive(model_ifstream);
      input_archive >> surfpackCModel; 
      std::cout << "Model loaded from binary file '" << model_filename << "'." 
		<< std::endl;
    }
    else {
      std::cout << "Loading model text file '" << model_filename << "'." 
		<< std::endl;
      boost::archive::text_iarchive input_archive(model_ifstream);
      input_archive >> surfpackCModel; 
      std::cout << "Model loaded from text file '" << model_filename << "'." 
		<< std::endl;
    }
 
  }
  catch (const std::exception& e) {
    std::cerr << "Error loading surfpack model! Exception:\n" << e.what() 
	      << std::endl;
    success = false;
  } 
  catch (const std::string& s) {
    std::cerr << "Error loading model surfpack! String:\n" << s << std::endl;
    success = false;
  }
  catch (...) {
    std::cerr << "Error loading surfpack model! Unknown error." << std::endl;
    success = false;
  }

  return (success ? 0 : 1);
}


extern "C"
double surfpack_eval_model(const double * const eval_pt, unsigned int num_vars)
{
  // evaluate the surfpack model, using a std::vector
  std::vector<double> eval_vec(&eval_pt[0], &eval_pt[num_vars-1]);
  return surfpackCModel->operator()(eval_vec);
}


extern "C"
void surfpack_free_model()
{
  delete surfpackCModel;
}

