#include "config.h"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>
#include <set>
#include "vector_enhancements.h"
#include "ann.h"
#include "surfpack.h"
#include "SurfPoint.h"
#include "SurfData.h"
#include "Surface.h"
#include "ANNSurface.h"

extern "C" void dgels_(char& trans, int& m, int& n, int& nrhs, double* A,
      int& lda, double* B, int& ldb, double* work, int& lwork, int& info);

extern "C" void dgemm_(char& transa, char& transb, int& m, int& n, int& k, 
  double& alpha, double* A, int& lda, double* B, int& ldb, double& beta, 
  double* C, int& ldc);

extern "C" void dgemv_(char& trans, int& m, int& n, double& alpha, 
  double* A, int& lda, double* x, int& incx, double& beta, double* y, 
  int& incy);

using namespace std;
//_____________________________________________________________________________
// Data members 
//_____________________________________________________________________________

const string ANNSurface::name = "ANN";

//_____________________________________________________________________________
// Creation, Destruction, Initialization
//_____________________________________________________________________________

ANNSurface::ANNSurface(SurfData* sd) : Surface(sd),
  annObject(0)
{
#ifdef __TESTING_MODE__
  constructCount++;
#endif
  init();
}

ANNSurface::ANNSurface(const string filename) : Surface(0),
  annObject(0)
{
#ifdef __TESTING_MODE__
  constructCount++;
#endif
  read(filename);
}
ANNSurface::~ANNSurface()
{
  delete annObject;
#ifdef __TESTING_MODE__
  destructCount++;
#endif
}

void ANNSurface::init()
{
  norm_bound = 0.8; 
  percent = 0; 
  svdfactor = 0.90;
}

//_____________________________________________________________________________
// Overloaded Operators 
//_____________________________________________________________________________

//_____________________________________________________________________________
// Queries
//_____________________________________________________________________________

const std::string ANNSurface::surfaceName() const
{
  return name;
}

unsigned ANNSurface::minPointsRequired() const
{
  if (xsize <= 0) {
    throw string(
      "Dimensionality of data needed to determine number of required samples."
    );
  } else {
    return 1+xsize;
  }
}

double ANNSurface::evaluate(const std::vector<double>& x)
{
  vector<double> annOutput(1);
  annObject->map(x, &annOutput);
  return annOutput[0];
}

//_____________________________________________________________________________
// Commands 
//_____________________________________________________________________________

void ANNSurface::build(SurfData& data)
{
  delete annObject;
  annObject = new ANNApprox;
  vector< vector< double > > training_inputs;
  vector< vector< double > > training_outputs;
  cout << "Points: " << data.size() << " Dimensions: " << data.xSize() << endl;
  reshape_2d(training_inputs, data.size(), data.xSize());
  reshape_2d(training_outputs, data.size(), 1);
  for (unsigned i = 0; i < data.size(); i++) {
    const SurfPoint& sp = data[i];
    for (unsigned j = 0; j < data.xSize(); j++) {
      training_inputs[i][j] = sp.X()[j];
      //cout << "inputs[" << i << "][" << j << "]: " << training_inputs[i][j] << endl;

      training_outputs[i][0] = data.getResponse(i);
      //cout << "outputs[" << i << "]: " << training_outputs[i][0] << endl;
    }
  }
  double local_norm_bound = 0.8, local_percent = 0, local_svdfactor = 0.99;
  annObject->normalize_data(training_inputs, training_outputs, local_norm_bound);
  annObject->set_aside_test_exemplars(local_percent);
  int num_neurons = annObject->numExemplars - 1;
  annObject->build_approximation(local_svdfactor, num_neurons);

  //for (unsigned k = 0; k < data.size(); k++) {
  //  cout << "Prediction " << k << " " << evaluate(data[k].X()) << endl;
  //}
}

void ANNSurface::config(const SurfpackParser::Arg& arg)
{
  string argname = arg.name;
  if (argname == "norm_bound") {
    norm_bound = arg.rval.real;
  } else if (argname == "svdfactor") {
    svdfactor = arg.rval.real;
  } else if (argname == "fraction_withheld") {
    percent = arg.rval.real;
  } else {
    Surface::config(arg);
  }
}

/// Create a surface of the same type as 'this.'  This objects data should
/// be replaced with the dataItr passed in, but all other attributes should
/// be the same (e.g., a second-order polynomial should return another 
/// second-order polynomial.  Surfaces returned by this method can be used
/// to compute the PRESS statistic.
ANNSurface* ANNSurface::makeSimilarWithNewData(SurfData* surfData)
{
  return new ANNSurface(surfData);
}

//_____________________________________________________________________________
// Helper methods 
//_____________________________________________________________________________

//_____________________________________________________________________________
// I/O 
//_____________________________________________________________________________

void ANNSurface::writeBinary(std::ostream& os)
{
  if (annObject) annObject->writeBinary(os);
}

void ANNSurface::writeText(std::ostream& os)
{
  if (annObject) annObject->writeText(os);

}

void ANNSurface::readBinary(std::istream& is)
{

  delete annObject;
  annObject = new ANNApprox;
  annObject->readBinary(is);
}

void ANNSurface::readText(std::istream& is)
{
  delete annObject;
  annObject = new ANNApprox;
  annObject->readText(is);
}

//_____________________________________________________________________________
// Testing 
//_____________________________________________________________________________

#ifdef __TESTING_MODE__
  int ANNSurface::constructCount = 0;
  int ANNSurface::destructCount = 0;
#endif
