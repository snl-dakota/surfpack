#include "SurfpackModel.h"
#include "SurfData.h"

using std::cout;
using std::endl;
using std::vector;


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
  : params(), ndims(0)
{

}

SurfpackModelFactory::SurfpackModelFactory(const ParamMap& params_in)
  : params(params_in), ndims(0)
{

}

void SurfpackModelFactory::config()
{
  ndims = atoi(params["ndims"].c_str());
  assert(ndims);
}

unsigned SurfpackModelFactory::minPointsRequired()
{
  SurfpackModelFactory::config();
  assert(ndims);
  return (ndims+1)*(ndims+2)/2; 
}

const ParamMap& SurfpackModelFactory::parameters() const
{
  return params;
}

void SurfpackModelFactory::add(std::string name, std::string value)
{
  params.insert(ModelParam(name,value));
}


