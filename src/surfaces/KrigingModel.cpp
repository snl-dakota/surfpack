#include "surfpack_system_headers.h"
#include "KrigingModel.h"
#include "SurfpackMatrix.h"
#include "SurfData.h"
#include "surfpack.h"
#include "ModelScaler.h"
#include "AxesBounds.h"
#include "ModelFitness.h"

using std::cout;
using std::endl;
using std::copy;
using std::vector;
using std::string;


void KrigingModel::
surfdata_to_nkm_surfdata(const SurfData& sd, nkm::SurfData& nkm_sd)
{
  unsigned num_points = sd.size();
  unsigned x_size = sd.xSize();
  unsigned f_size = sd.fSize();

  nkm::MtxDbl XR(num_points, x_size), Y(num_points, f_size);

  for (unsigned point_index=0; point_index<num_points; ++point_index) {
    const SurfPoint& point = sd[point_index];
    const VecDbl x = point.X();
    for (unsigned x_index=0; x_index<x_size; ++x_index)
      XR(point_index, x_index) = x[x_index];
    for (unsigned f_index=0; f_index<f_size; ++f_index)
      Y(point_index, f_index) = point.F(f_index);
  }

  nkm_sd = nkm::SurfData(XR, Y);
}


KrigingModel::KrigingModel(const SurfData& sd, const ParamMap& args)
  : SurfpackModel(sd.xSize())
{
  surfdata_to_nkm_surfdata(sd, nkmSurfData);
  nkmKrigingModel = new nkm::KrigingModel(nkmSurfData, args);
  nkmKrigingModel->create();
}


KrigingModel::~KrigingModel()
{
  delete nkmKrigingModel;
}


double KrigingModel::evaluate(const VecDbl& x) const
{
  nkm::MtxDbl nkm_x(1, ndims);
  for(size_t i=0; i<ndims; ++i)
    nkm_x(0, i) = x[i];
  return (nkmKrigingModel->evaluate(nkm_x));
}


double KrigingModel::variance(const VecDbl& x) const
{
  nkm::MtxDbl nkm_x(1, ndims);
  for(size_t i=0; i<ndims; ++i)
    nkm_x(0, i) = x[i];
  //double adj_var=nkmKrigingModel->eval_variance(nkm_x);
  //printf("KrigingModel::variance() adj_var=%g\n",adj_var);
  //return (adj_var);
  return (nkmKrigingModel->eval_variance(nkm_x));
}


VecDbl KrigingModel::gradient(const VecDbl& x) const
{
  nkm::MtxDbl nkm_x(1, ndims);
  for(int i=0; i<ndims; ++i)
    nkm_x(0, i) = x[i];

  nkm::MtxDbl nkm_d1y(1, ndims);
  nkmKrigingModel->evaluate_d1y(nkm_d1y, nkm_x);

  VecDbl d1y(ndims, 0.0); 
  for(int i=0; i<ndims; ++i)
    d1y[i] = nkm_d1y(i);

  return d1y;
}

MtxDbl KrigingModel::hessian(const VecDbl& x) const
{
  nkm::MtxDbl nkm_x(1, ndims);
  for(int i=0; i<ndims; ++i)
    nkm_x(0, i) = x[i];

  int num_lower_elem=(ndims+1)*ndims/2;
  nkm::MtxDbl nkm_d2y(1, num_lower_elem);
  nkmKrigingModel->evaluate_d2y(nkm_d2y, nkm_x);

  MtxDbl d2y(ndims, ndims, 0.0); 
  int k=0;
  for(int j=0; j<ndims; ++j) {
    d2y(j,j)=nkm_d2y(k);
    k++;
    for(int i=j+1; i<ndims; ++i) {
      d2y(i,j) = nkm_d2y(k);
      d2y(j,i) = d2y(i,j);
      k++;
    }
  }

  return d2y;
}


std::string KrigingModel::asString() const
{
  return (nkmKrigingModel->model_summary_string());
  // TODO: be able to write NKM_KrigingModel as a string
  //assert(false);
}


//KrigingModel KrigingModel::Create(const SurfData& sd)
//{
//  ConminKriging ck(sd);
//  VecDbl best_guess(sd.xSize(),1.0);
//  unsigned max_iter = 100; 
//  double opt_val;
//  VecDbl rhs;
//  ck.optimize(best_guess,opt_val,max_iter); 
//  KrigingBasisSet kbs(SurfData::asVecVecDbl(sd),best_guess);
//}

///////////////////////////////////////////////////////////
/// 	Kriging Model Factory	
///////////////////////////////////////////////////////////

KrigingModelFactory::KrigingModelFactory()
  : SurfpackModelFactory()
{

}

KrigingModelFactory::KrigingModelFactory(const ParamMap& args)
  : SurfpackModelFactory(args)
{

}

void KrigingModelFactory::config()
{
  SurfpackModelFactory::config();
}

SurfpackModel* KrigingModelFactory::Create(const std::string& model_string)
{
  ///\todo Be able to parse a Kriging model from a string
  assert(false);
  return 0;
}


typedef std::pair<double,VecDbl> KMPair;
SurfpackModel* KrigingModelFactory::Create(const SurfData& sd)
{
  this->add("ndims",surfpack::toString(sd.xSize()));
  this->config();

  return new KrigingModel(sd, params);
}
