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
  std::vector<std::vector<nkm::MtxDbl> > derY(f_size);
  nkm::MtxInt der_order(1,f_size); 
  der_order.zero(); //set contents to zero

  if(num_points>0) {
    //could increment der_order independently for each dimension but old surfdata/surfpoint does not support this, nkm::SurfData supports arbitrarily high order derivatives.

    const SurfPoint& point=sd[0];
    if(point.fGradientsSize() > 0) {
      for(unsigned f_index=0; f_index<f_size; ++f_index)
	++der_order(f_index);
      if(point.fHessiansSize() > 0)
	for(unsigned f_index=0; f_index<f_size; ++f_index)
	  ++der_order(f_index);
    }
    
    for(unsigned f_index=0; f_index<f_size; ++f_index) {
      derY[f_index].resize(der_order(f_index)+1);
      for(unsigned der_order_index=1; der_order_index<=der_order(f_index); ++der_order_index)
	derY[f_index][der_order_index].newSize(num_points,nkm::num_multi_dim_poly_coef(x_size,-der_order_index));
    }
  }


  for (unsigned point_index=0; point_index<num_points; ++point_index) {
    const SurfPoint& point = sd[point_index];
    const VecDbl x = point.X();
    for (unsigned x_index=0; x_index<x_size; ++x_index)
      XR(point_index, x_index) = x[x_index];
    for (unsigned f_index=0; f_index<f_size; ++f_index) 
      Y(point_index, f_index) = point.F(f_index);


    // given NKM ordering of derY, probably need to loop differently
    // to populate the matrices; this just for demo

    // example of mapping first derivatives
    // there should be 0 or f_size gradients (could throw error)
    if (point.fGradientsSize() > 0) {
      for (unsigned f_index=0; f_index < f_size; ++f_index) {
	assert(der_order(f_index)>=1);  //could change this to a throw
	const vector<double>& sd_gradient = point.fGradient(f_index);
	//cout << "Surfpack gradient for point " << point_index << ", function "
	//     << f_index << ": [ ";
	for (unsigned x_index=0; x_index < x_size; ++x_index) {
	  // accessing each gradient element
	  //cout << sd_gradient[x_index] << " ";
	  derY[f_index][1](point_index,x_index) =sd_gradient[x_index];
	}
	//cout << "]" << endl;
      }
    }
    else{
      for (unsigned f_index=0; f_index < f_size; ++f_index) 
	assert(der_order(f_index)==0);  //could change this to a throw
    }

    // example of mapping second derivatives
    // there should be 0 or f_size Hessians (could throw error)
    if (point.fHessiansSize() > 0) {
      for (unsigned f_index=0; f_index<f_size; ++f_index) {
	assert(der_order(f_index)>=2);  //could change this to a throw

	const SurfpackMatrix<double>& sd_hessian = point.fHessian(f_index);
	//cout << "Surfpack Hessian for point " << point_index << ", function "
	//     << f_index << " is:\n";
	unsigned der_index=0;
	for (unsigned xj_index=0; xj_index < x_size; ++xj_index) 
	  for (unsigned xk_index=xj_index; xk_index < x_size; ++xk_index, ++der_index)
	    derY[f_index][2](point_index,der_index)=
	      sd_hessian(xj_index,xk_index);

	//for (unsigned xj_index=0; xj_index < x_size; ++xj_index) {
	//  for (unsigned xk_index=0; xk_index < x_size; ++xk_index) {
	//    // accessing each Hessian element
	//    cout << sd_hessian(xj_index, xk_index) << " ";
	//  }
	//  cout << "\n";
	//}
	//cout << endl;
      }
    }
    else{
      for(unsigned f_index=0; f_index<f_size; ++f_index)
	assert(der_order(f_index)<2);
    }
  }

  // TODO: populate with derY as well
  nkm_sd = nkm::SurfData(XR, Y, der_order, derY);
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
