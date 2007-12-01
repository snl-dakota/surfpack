#include "surfpack_system_headers.h"
#include "surfpack.h"
#include "LinearRegressionModel.h"
#include "SurfData.h"


using std::cout;
using std::endl;
using std::string;

double LRMBasisSet::eval(unsigned index, const VecDbl& x) const
{
  assert(index < bases.size());
  double result = 1.0;
  for(VecUnsIt it = bases[index].begin(); it != bases[index].end(); ++it) {
    if (*it >= x.size()) cout << index << " " << *it << endl;
    assert(*it < x.size());
    result *= x[*it];
  } 
  return result;
}

double LRMBasisSet::deriv(unsigned index, const VecDbl& x, const VecUns& vars) const
{
  VecUns counts(x.size(),0);
  for(VecUnsIt it = bases[index].begin(); it != bases[index].end(); ++it) {
    assert(*it < x.size());
    counts[*it]++;
  } 
  double coeff = 1.0;
  for(VecUns::const_iterator it = vars.begin(); it != vars.end(); ++it) {
    assert(*it < x.size());
    unsigned index = *it;
    // Taken derivative with respect to this variable too many times
    if (!counts[index]) return 0.0; 
    coeff *= counts[index]--;
  }
  unsigned sum = accumulate(counts.begin(),counts.end(),0);
  if (sum == 0) return coeff; // no vars left
  double term = 1.0;
  for(unsigned v = 0; v < counts.size(); v++) {
    for (unsigned c = 0; c < counts[v]; c++) {
      term *= x[v];
    }
  }
  //cout << "Deriv: " << index << " X: ";
  //copy(x.begin(),x.end(),std::ostream_iterator<double>(cout," "));
  //cout << " Vars: " ;
  //copy(vars.begin(),vars.end(),std::ostream_iterator<double>(cout," "));
  //cout << " Val: " << coeff*term << endl;
  return coeff*term;
}
  
std::string LRMBasisSet::asString() const
{
  std::ostringstream os;
  for(VecVecUns::const_iterator it = bases.begin(); it != bases.end(); ++it) {
    if (it->empty()) { os << "Unity\n"; continue; }
    copy(it->begin(),it->end(),std::ostream_iterator<unsigned>(os," "));
    os << "\n";
  }
  return os.str();
}

void LRMBasisSet::add(const std::string& s_basis)
{
  bases.push_back(surfpack::toVec<unsigned>(s_basis));
}

LinearRegressionModel::LinearRegressionModel(const unsigned dims, 
  const LRMBasisSet& bs_in, const VecDbl& coeffs_in)
  : SurfpackModel(dims), bs(bs_in), coeffs(coeffs_in)
{
  assert(bs.bases.size() == coeffs.size());
}

double LinearRegressionModel::evaluate(const VecDbl& x) const
{
  assert(coeffs.size() == bs.bases.size());
  double sum = 0;
  for (unsigned i = 0; i < coeffs.size(); i++) {
    sum += coeffs[i]*bs.eval(i,x);
  }
  return sum;
}

VecDbl LinearRegressionModel::gradient(const VecDbl& x) const
{
  assert(!x.empty());
  cout << "IN gradient x[0] = " << x[0] << endl;
  assert(coeffs.size() == bs.bases.size());
  VecUns diff_var(1,0); // variable with which to differentiate
  VecDbl result(x.size(),0.0);
  for (unsigned i = 0; i < x.size(); i++) {
    diff_var[0] = i;
    for (unsigned j = 0; j < bs.bases.size(); j++) {
      result[i] += coeffs[j]*bs.deriv(j,x,diff_var);
    }
  }
  return result;
}

std::string LinearRegressionModel::asString() const
{
  std::ostringstream os;
  os << "\nbases:\n" << bs.asString() << "coeffs: ";
  copy(coeffs.begin(),coeffs.end(),std::ostream_iterator<double>(os," "));
  os << "\n";
  return os.str();
}

///////////////////////////////////////////////////////////
///	Moving Least Squares Model Factory
///////////////////////////////////////////////////////////

VecDbl LinearRegressionModelFactory::lrmSolve(const LRMBasisSet& bs, const ScaledSurfData& ssd)
{
  MtxDbl A(ssd.size(),bs.size(),true);
  for (unsigned i = 0; i < ssd.size(); i++) {
    for (unsigned j = 0; j < bs.size(); j++) {
      A(i,j) = bs.eval(j,ssd(i));
    }
  }
  VecDbl b = ssd.getResponses();
  VecDbl x;
  if (eqConRHS.empty()) {
    surfpack::linearSystemLeastSquares(A,x,b);
  } else {
    surfpack::leastSquaresWithEqualityConstraints(A,x,b,eqConLHS,eqConRHS);
  }
  return x; 
}

LRMBasisSet LinearRegressionModelFactory::CreateLRM(unsigned order, 
  unsigned dims)
{
  LRMBasisSet bs;
  bs.add(std::string(""));
  //unsigned order = 3;
  //unsigned dims = 3;
  std::deque<Term> q;
  q.push_front(Term(VecUns()));
  while (!q.empty()) {
    Term& t = q.front();
    VecUns& v = t.vars;
    if (v.size() < order && !t.color) { // extendable
      t.color = true; // only extend it once
      Term new_term = Term(v);
      if (v.empty()) new_term.vars.push_back(0);
      else new_term.vars.push_back(v.back());
      //std::copy(new_term.vars.begin(),new_term.vars.end(),std::ostream_iterator<unsigned>(std::cout," ")); std::cout << "\n";
      bs.bases.push_back(new_term.vars);
      q.push_front(new_term);
    } else if (!v.empty() && v.back() < dims-1) {
      v.back()++;
      //std::copy(v.begin(),v.end(),std::ostream_iterator<unsigned>(std::cout," ")); std::cout << "\n";
      bs.bases.push_back(v);
      t.color = false;
    } else {
      q.pop_front();
    }
  }
  return bs;
} 

SurfpackModel* LinearRegressionModelFactory::Create(const std::string& model_string)
{
  ///\todo Be able to parse an RBF model from a string
  assert(false);
  return 0;
}

SurfpackModel* LinearRegressionModelFactory::Create(const SurfData& sd)
{
  this->add("ndims",surfpack::toString(sd.xSize()));
  this->config();
  ModelScaler* ms = NormalizingScaler::Create(sd);
  ScaledSurfData ssd(*ms,sd);
  
  LRMBasisSet bs = CreateLRM(2,sd.xSize());
  cout << bs.asString() << endl;
  cout << "sd size: " << sd.size() << endl;
  VecDbl coeffs = lrmSolve(bs,ssd);
  copy(coeffs.begin(),coeffs.end(),std::ostream_iterator<double>(cout,"|"));
  cout << "\n";
  SurfpackModel* lrm = new LinearRegressionModel(sd.xSize(),bs,coeffs);
  lrm->scaler(ms);
  delete ms;
  return lrm;
}

unsigned LinearRegressionModelFactory::minPointsRequired()
{
  config();
  LRMBasisSet bs = CreateLRM(order,ndims);
  return bs.size() - eqConRHS.size();
}

LinearRegressionModelFactory::LinearRegressionModelFactory()
  : SurfpackModelFactory(), order(2)
{

}

LinearRegressionModelFactory::LinearRegressionModelFactory(const ParamMap& args)
  : SurfpackModelFactory(args), order(2)
{

}

void LinearRegressionModelFactory::config()
{
  SurfpackModelFactory::config();
  string strarg;
  strarg = params["order"];
  if (strarg != "") order = atoi(strarg.c_str());
}

void LinearRegressionModelFactory::setEqualityConstraints(unsigned asv,const SurfPoint& sp,  double valuePtr, VecDbl& gradient, MtxDbl& hessian)
{
  config();
  LRMBasisSet bs = CreateLRM(order,ndims);
  VecDbl coefficients(bs.size());
  unsigned numConstraints = 0;
  if (asv & 1) numConstraints += 1; // value at a particular point
  if (asv & 2) numConstraints += ndims; // gradient at a point
  if (asv & 4) numConstraints += (ndims*ndims+ndims)/2; // hessian at a point
  eqConRHS.resize( numConstraints );
  // Must compute number of terms first
  //MtxDbl temp(eqConRHS.size(),coefficients.size(),true);
  //eqConLHS = temp;
  eqConLHS.reshape(eqConRHS.size(),coefficients.size());
  // Marks the index of the next constraint to be added (necessary since
  // indices of e.g. the gradient constraints will be different depending on
  // whether or not the value constraint is used
  unsigned index = 0;
  // If requested, add the equality constraint for the point value
  if (asv & 1) {
    for (unsigned i = 0; i < bs.size(); i++) {
      eqConLHS(index,i) = bs.eval(i,sp.X());
      eqConRHS[index] = valuePtr;
      ++index;
    }
  }

  // If requested, add the equality constraints for the gradient
  if (asv & 2) {
    //const VecDbl& gradient = *gradientPtr;
    assert(gradient.size() == ndims);
    VecUns factorCounts;
    VecUns diff_counts;
    for (unsigned dif_var = 0; dif_var < ndims; dif_var++ ) {
      diff_counts = VecUns(ndims,0);
      diff_counts[dif_var] = 1;
      for (unsigned i = 0; i < bs.size(); i++) {
        eqConLHS(index,i) = bs.deriv(i,sp.X(), diff_counts);
        ++index;
      }
      eqConRHS[index] = gradient[dif_var];
      ++index;
    }
  }

  // If requested, add the equality constraints for the hessian
  if (asv & 4) {
    //MtxDbl& hessian = *hessianPtr;
    assert(hessian.getNCols() == ndims);
    assert(hessian.getNRows() == ndims);
    VecUns factorCounts;
    VecUns diff_counts;
    for (unsigned dif_var1 = 0; dif_var1 < ndims; dif_var1++ ) {
      for (unsigned dif_var2 = dif_var1; dif_var2 < ndims; dif_var2++ ) {
        diff_counts = VecUns(ndims,0);
        diff_counts[dif_var1]++;
        diff_counts[dif_var2]++;
        for (unsigned i = 0; i < bs.size(); i++) {
          eqConLHS(index,i) = bs.deriv(i,sp.X(), diff_counts);
          ++index;
        }
        eqConRHS[index] = hessian(dif_var1,dif_var2);
        ++index;
      } // dif_var2
    } // dif_var1
  } // if hessian needed
}
