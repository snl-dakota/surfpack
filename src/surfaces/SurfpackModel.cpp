#include "SurfpackModel.h"
#include "SurfData.h"
#include "surfpack.h"

using std::cout;
using std::endl;
using std::vector;
using std::string;


///////////////////////////////////////////////////////////
///	Surfpack Model 
///////////////////////////////////////////////////////////

VecDbl SurfpackModel::operator()(const SurfData& data) const
{
  VecDbl result(data.size());
  for (unsigned pt = 0; pt < data.size(); pt++) {
    result[pt] = (*this)(data(pt));
  }
  return result;
}

double SurfpackModel::operator()(const VecDbl& x) const
{
  //cout << "\nunscaled\n";
  //copy(x.begin(),x.end(),std::ostream_iterator<double>(cout," "));
  const VecDbl& x1 = mScaler->scale(x);
  //cout << "\nscaled\n";
  //copy(x1.begin(),x1.end(),std::ostream_iterator<double>(cout," "));
  double value = evaluate(x1);
  //cout << "evaluated: " << value << "\n";
  double result = mScaler->descale(value);
  //cout << "descaled: " << result << endl;
  return result;
  
  //return mScaler->descale(evaluate(mScaler->scale(x)));
}

SurfpackModel::SurfpackModel(unsigned ndims_in) 
  : ndims(ndims_in), mScaler(new NonScaler)
{

}

SurfpackModel::SurfpackModel(const SurfpackModel& other)
  : ndims(other.ndims), mScaler(other.mScaler->clone())
{

}

SurfpackModel::~SurfpackModel()
{
  delete mScaler; mScaler = 0;
}

VecDbl SurfpackModel::gradient(const VecDbl& x) const
{
  throw std::string("This model does not currently support gradients");
}

MtxDbl SurfpackModel::hessian(const VecDbl& x) const
{
  throw std::string("This model does not currently support hessians");
}

void SurfpackModel::scaler(ModelScaler* ms)
{
  delete mScaler;
  mScaler = ms->clone();
}

ModelScaler* SurfpackModel::scaler() const
{
  return mScaler;
}

///////////////////////////////////////////////////////////
///	Surfpack Model Factory
///////////////////////////////////////////////////////////

SurfpackModelFactory::SurfpackModelFactory()
  : params(), ndims(0), response_index(0)
{

}

SurfpackModelFactory::SurfpackModelFactory(const ParamMap& params_in)
  : params(params_in), ndims(0), response_index(0)
{

}

void SurfpackModelFactory::config()
{
  ndims = atoi(params["ndims"].c_str());
  assert(ndims);
  string arg = params["response_index"];
  if (arg != "") response_index = atoi(arg.c_str());
}

unsigned SurfpackModelFactory::minPointsRequired()
{
  SurfpackModelFactory::config();
  assert(ndims);
  return (ndims+1);
}

unsigned SurfpackModelFactory::recommendedNumPoints()
{
  SurfpackModelFactory::config();
  assert(ndims);
  return (5*ndims);
}

const ParamMap& SurfpackModelFactory::parameters() const
{
  return params;
}

void SurfpackModelFactory::add(const std::string& name, const std::string& value)
{
  params.insert(ModelParam(name,value));
}

SurfpackModel* SurfpackModelFactory::Build(const SurfData& sd)
{
  //cout << "Data:\n";
  //sd.writeText(cout);
  //cout << "\nParams:\n";
  //
  //  for (ParamMap::iterator itr = params.begin();
  //      itr != params.end(); itr++) {
  //    std::cout << "     " << itr->first << ": " << itr->second << std::endl;
  //  }
  this->add("ndims",surfpack::toString<unsigned>(sd.xSize()));
  this->config();
  std::cout << "SurfpackModelFactory built with parameters:" << std::endl; 
  for (ParamMap::iterator itr = params.begin();
       itr != params.end(); itr++) {
    std::cout << "     " << itr->first << ": " << itr->second << std::endl;
  }
  sd.setDefaultIndex(this->response_index);
  if (sd.size() < minPointsRequired()) {
    throw string("Not enough Points");
  }
  SurfpackModel* model = Create(sd);
  model->parameters(params);
  return model;
}
