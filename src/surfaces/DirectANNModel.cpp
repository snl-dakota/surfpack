#include "surfpack_system_headers.h"
#include "DirectANNModel.h"
#include "SurfData.h"
#include "surfpack.h"


using std::cout;
using std::endl;
using std::string;

DirectANNBasisSet::DirectANNBasisSet(const MtxDbl& weights_in)
  : weights(weights_in)
{

}

double DirectANNBasisSet::nodeSum(unsigned index, const VecDbl& x) const
{
  assert(index < weights.getNRows());
  assert(x.size() + 1 == weights.getNCols());
  double sum = 0.0;
  for (unsigned i = 0; i < x.size(); i++) {
    sum += weights(index,i)*x[i];
  }
  sum += weights(index,x.size()); // bias weight
  return sum;
}

double DirectANNBasisSet::eval(unsigned index, const VecDbl& x) const
{
  //printf("sum: %f tanh thereof: %f\n",nodeSum(index,x),tanh(nodeSum(index,x)));
  return tanh(nodeSum(index,x));
}

double DirectANNBasisSet::deriv(unsigned index, const VecDbl& x, const VecUns& vars) const
{
  assert(vars.size() == 1);
  assert(vars[0] < x.size());
  double sum = nodeSum(index,x);
  double tanhsum = tanh(sum);
  return (1.0 - tanhsum*tanhsum)*weights(index,vars[0]);
}
  
std::string DirectANNBasisSet::asString() const
{
  return weights.asString();
}

DirectANNModel::DirectANNModel(const DirectANNBasisSet& bs_in, const VecDbl& coeffs_in)
  : SurfpackModel(bs_in.weights.getNCols()), bs(bs_in), coeffs(coeffs_in)
{
  assert(bs.weights.getNRows()+1 == coeffs.size());
}

double DirectANNModel::evaluate(const VecDbl& x) const
{
  assert(coeffs.size() == bs.weights.getNRows() + 1);
  double sum = 0;
  //cout << "-----> x: ";
  //copy(x.begin(),x.end(),std::ostream_iterator<double>(cout," "));
  //cout << "\n";
  for (unsigned i = 0; i < bs.weights.getNRows(); i++) {
    //printf("  i: %d coeff: %f node: %f\n",i,coeffs[i],bs.eval(i,x));
    sum += coeffs[i]*bs.eval(i,x);
  }
  sum += coeffs.back(); // bias weight 
  //printf("  sum: %f bias: %f tanh: %f\n",sum,coeffs.back(),tanh(sum));
  return tanh(sum);
}

VecDbl DirectANNModel::gradient(const VecDbl& x) const
{
  assert(!x.empty());
  assert(x.size() + 1 == bs.weights.getNCols());
  VecUns diff_var(1,0); // variable with which to differentiate
  VecDbl nodeSums(bs.weights.getNRows());
  double finalSum; // the unsigmoided value of the output node
  for (unsigned r = 0; r < bs.weights.getNRows(); r++) {
    nodeSums[r] = bs.nodeSum(r,x);
    finalSum += coeffs[r]*tanh(nodeSums[r]);
  }
  double tanhsum = tanh(finalSum+coeffs[bs.weights.getNRows()]);
  double finalSumMultiplier = 1 - tanhsum*tanhsum;
  VecDbl result(x.size(),0.0);
  for (unsigned v = 0; v < x.size(); v++) {
    for (unsigned i = 0; i < bs.weights.getNRows(); i++) {
      double tanhNodeSum = tanh(nodeSums[i]);
      result[v] += coeffs[i]*(1-tanhNodeSum*tanhNodeSum)*bs.weights(i,v);
    }
    result[v] *= finalSumMultiplier;
  }
  return result;
}

std::string DirectANNModel::asString() const
{
  std::ostringstream os;
  os << "\nweights:\n" << bs.asString() << "coeffs: ";
  copy(coeffs.begin(),coeffs.end(),std::ostream_iterator<double>(os," "));
  os << "\n";
  return os.str();
}

MtxDbl randomMatrix(unsigned nrows, unsigned ncols)
{
  MtxDbl rm(nrows,ncols);
  for (unsigned i = 0; i < nrows; i++) {
    for (unsigned j = 0; j < ncols; j++) {
      rm(i,j) = (double)rand()/INT_MAX;
    }
  }
  return rm;
}

///////////////////////////////////////////////////////////
/// 	DirectANN Model Factory	
///////////////////////////////////////////////////////////

DirectANNModelFactory::DirectANNModelFactory()
  : SurfpackModelFactory(), nNodes(0)
{

}

DirectANNModelFactory::DirectANNModelFactory(const ParamMap& args)
  : SurfpackModelFactory(args), nNodes(0)
{

}

void DirectANNModelFactory::config()
{
  SurfpackModelFactory::config();
  string strarg;
  strarg = params["nodes"];
  if (strarg != "") nNodes = atoi(strarg.c_str()); 
}

SurfpackModel* DirectANNModelFactory::Create(const std::string& model_string)
{
  ///\todo Be able to parse an RBF model from a string
  assert(false);
  return 0;
}


typedef std::pair<double,VecDbl> KMPair;
SurfpackModel* DirectANNModelFactory::Create(const SurfData& sd)
{
  this->add("ndims",surfpack::toString(sd.xSize()));
  this->config();
  const unsigned maxnodes = 100;
  assert(sd.size());
  assert(sd.xSize());
  if (!nNodes) nNodes = sd.size();
  if (nNodes > maxnodes) nNodes = maxnodes;
  MtxDbl random_weights = randomMatrix(nNodes,sd.xSize()+1);
  DirectANNBasisSet bs(random_weights);
  MtxDbl A(sd.size(),nNodes+1,true);
  VecDbl b(sd.size(),0.0);
  for (unsigned samp = 0; samp < sd.size(); samp++) {
    for (unsigned n = 0; n < nNodes; n++) { 
      A(samp,n) = bs.eval(n,sd(samp));
    }
    A(samp,nNodes) = 1.0; // for hidden layer bias
    b[samp] = atanh(sd.getResponse(samp));
  }
  VecDbl x;
  surfpack::linearSystemLeastSquares(A,x,b);
  return new DirectANNModel(bs,x);
}
