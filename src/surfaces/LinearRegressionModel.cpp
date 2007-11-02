#include "surfpack_system_headers.h"
#include "surfpack.h"
#include "LinearRegressionModel.h"
#include "SurfData.h"

using std::cout;
using std::endl;

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

///////////////////////////////////////////////////////////////////////////////
//////////////           FACTORY METHODS     //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

VecDbl LinearRegressionModel::lrmSolve(const LRMBasisSet& bs, const ScaledSurfData& ssd)
{
  MtxDbl A(ssd.size(),bs.size(),true);
  for (unsigned i = 0; i < ssd.size(); i++) {
    for (unsigned j = 0; j < bs.size(); j++) {
      A(i,j) = bs.eval(j,ssd(i));
    }
  }
  VecDbl b = ssd.getResponses();
  VecDbl x;
  surfpack::linearSystemLeastSquares(A,x,b);
  return x; 
}

//LinearRegressionModel
LRMBasisSet LinearRegressionModel::CreateLRM(unsigned order, unsigned dims)
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

LinearRegressionModel LinearRegressionModel::Create(const SurfData& sd)
{
  ModelScaler* ms = NormalizingScaler::Create(sd);
  ScaledSurfData ssd(*ms,sd);
  
  LRMBasisSet bs = CreateLRM(2,sd.xSize());
  cout << bs.asString() << endl;
  cout << "sd size: " << sd.size() << endl;
  VecDbl coeffs = lrmSolve(bs,ssd);
  copy(coeffs.begin(),coeffs.end(),std::ostream_iterator<double>(cout,"|"));
  cout << "\n";
  LinearRegressionModel lrm(sd.xSize(),bs,coeffs);
  lrm.scaler(ms);
  delete ms;
  return lrm;
}
