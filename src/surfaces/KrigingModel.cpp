#include "surfpack_system_headers.h"
#include "KrigingModel.h"
#include "SurfpackMatrix.h"
#include "SurfData.h"
#include "surfpack.h"
#include "ModelScaler.h"
#include "AxesBounds.h"
#include "ModelFitness.h"
#include "Conmin.h"

using std::cout;
using std::endl;
using std::copy;
using std::vector;
using std::string;

double correlation_function(const VecDbl& correlations,
  const VecDbl& pt1, const VecDbl& pt2)
{
  double sum = 0.0;
  double mult_term = 0.0;
  for (unsigned i = 0; i < correlations.size(); i++) {
    mult_term = pt1[i]-pt2[i];
    sum += correlations[i]*mult_term*mult_term;
  }
  //cout << "sum: " << sum << " exp(-sum): " << exp(-sum) << endl;
  return exp(-sum);
}

KrigingBasisSet::KrigingBasisSet(const VecVecDbl& centers_in, const VecDbl& correlations_in)
  : centers(centers_in), correlations(correlations_in)
{

}

double KrigingBasisSet::eval(unsigned index, const VecDbl& x) const
{
  assert(!centers.empty());
  assert(x.size() == correlations.size());
  assert(x.size() == centers[0].size());
  double sum = 0.0;
  double diff = 0.0;
  for (unsigned i = 0; i < x.size(); i++) {
    diff = (x[i] - centers[index][i]);
    sum += correlations[i]*diff*diff;
  }
  return exp(-sum);
}

double KrigingBasisSet::deriv(unsigned index, const VecDbl& x, const VecUns& vars) const
{
  assert(vars.size() == 1);
  unsigned i = vars[0];
  return -2.0*eval(index,x)*correlations[i]*(x[i] - centers[index][i]);
}
  
std::string KrigingBasisSet::asString() const
{
  std::ostringstream os;
  os << "correlations: ";
  copy(correlations.begin(),correlations.end(),
       std::ostream_iterator<double>(os," "));
  os << "\npts:\n";
  for(VecVecDbl::const_iterator it = centers.begin(); it != centers.end(); ++it) {
    copy(it->begin(),it->end(),std::ostream_iterator<double>(os," "));
    os << "\n";
  }
  return os.str();
}

KrigingModel::KrigingModel(const KrigingBasisSet& bs_in, const VecDbl& rhs_in)
  : SurfpackModel(bs_in.correlations.size()), bs(bs_in), rhs(rhs_in)
{
  assert(bs.centers.size() == rhs.size());
}

KrigingModel::KrigingModel(const SurfData& sd, const VecDbl& correlations)
  : SurfpackModel(sd.xSize()), bs(SurfData::asVecVecDbl(sd),correlations)
{
  // Scale the data (in hopes of improving numerical properties)
  ModelScaler* ms = NonScaler::Create(sd);
  ScaledSurfData ssd(*ms,sd);
  this->scaler(ms);
  bs = KrigingBasisSet(ScaledSurfData::asVecVecDbl(ssd),correlations);

  MtxDbl R = getMatrix(ssd,correlations);
  //cout << "R(1,0): " << R.asString() << endl;
  VecDbl y = ssd.getResponses();
  
  // Invert the correlation matrix
  vector<int> ipvt;
  surfpack::LUFact(R,ipvt);
  //cout << R.asString() << endl;
  // Calculate the determinant before doing the inversion
  double determinantCorrMatrix = 1.0;
  for (unsigned i = 0; i < ssd.size(); i++) {
    determinantCorrMatrix *= R(i,i);
  }
  // perform the inverse operation in place and then rename the var
  surfpack::inverseAfterLUFact(R,ipvt);
  MtxDbl& Rinv = R; // just renaming the variable so as to not get confused
  vector<double> Rinv_times_y;
  surfpack::matrixVectorMult(Rinv_times_y,Rinv,y);
  double denominator_sum = Rinv.sum();
  double numerator_sum = surfpack::sum_vector(Rinv_times_y);
  this->betaHat = numerator_sum / denominator_sum;
  //cout << "betaHat: " << betaHat << endl;
  //cout << "numerator: " << numerator_sum << endl;
  //cout << "denominator: " << denominator_sum << endl;
  surfpack::vectorShift(y,betaHat);
  VecDbl& y_minus_betahat = y;
  surfpack::matrixVectorMult(this->rhs,Rinv,y_minus_betahat);
  double estVariance = surfpack::dot_product(y_minus_betahat,rhs);
  //cout << "EstVariance: " << estVariance << endl;
  //cout << "determinantCorrMatrix: " << determinantCorrMatrix << endl;
  //double likelihood = 
  //  -0.5*sd.size()*log(estVariance)+log(determinantCorrMatrix);
  this->likelihood = log(estVariance)+log(determinantCorrMatrix);
  cout << "Our likelihood: " << likelihood << endl;
  delete ms;
  
}

MtxDbl KrigingModel::getMatrix(const ScaledSurfData& ssd, const VecDbl& correlations) 
{
  unsigned n = ssd.size();
  MtxDbl R(n,n,true); // correlation matrix
  for (unsigned i = 0; i < n; i++) {
    for (unsigned j = 0; j < n; j++) {
      if (i == j) {
        R(i,j) = 1.0;
      } else {
        const VecDbl pt1 = ssd(i);
        const VecDbl& pt2 = ssd(j);
        R(i,j) = 
   	  correlation_function(correlations,pt1,pt2);
      }
    }
  }
  return R;
}

MtxDbl KrigingModel::corrMtx(const VecDbl& corr_vec, const SurfData& data)
{
  MtxDbl R(data.size(),data.size(),true); // correlation matrix
  for (unsigned i = 0; i < data.size(); i++) {
    for (unsigned j = 0; j < data.size(); j++) {
      if (i == j) {
        R(i,j) = 1.0;
      } else {
        R(i,j) = 
   	  correlation_function(corr_vec,data[i].X(),data[j].X());
      }
    }
  }
  return R;
}

double KrigingModel::evaluate(const VecDbl& x) const
{
  assert(rhs.size() == bs.centers.size());
  double sum = 0;
  for (unsigned i = 0; i < rhs.size(); i++) {
    sum += rhs[i]*bs.eval(i,x);
  }
  return betaHat+sum;
}

VecDbl KrigingModel::gradient(const VecDbl& x) const
{
  assert(!x.empty());
  cout << "IN gradient x[0] = " << x[0] << endl;
  assert(rhs.size() == bs.centers.size());
  VecUns diff_var(1,0); // variable with which to differentiate
  VecDbl result(x.size(),0.0);
  for (unsigned i = 0; i < x.size(); i++) {
    diff_var[0] = i;
    for (unsigned j = 0; j < bs.centers.size(); j++) {
      result[i] += rhs[j]*bs.deriv(j,x,diff_var);
    }
  }
  return result;
}

std::string KrigingModel::asString() const
{
  std::ostringstream os;
  os << "\ncenters:\n" << bs.asString() << "rhs: ";
  copy(rhs.begin(),rhs.end(),std::ostream_iterator<double>(os," "));
  os << "\n";
  return os.str();
}


ConminKriging::ConminKriging(const SurfData& data_in)
  : Conmin(data_in.xSize()), data(data_in)
{
  
}

void ConminKriging::optimize(VecDbl& x, double& final_val, unsigned max_iter)
{

  /// Size variable for CONMIN arrays. See CONMIN manual.
  /** N1 = number of variables + 2 */
  //int N1 = xsize + 2;
  ///// Size variable for CONMIN arrays. See CONMIN manual.
  ///** N2 = number of constraints + 2*(number of variables) */
  //int N2 = xsize * 2; // + NCON;
  ///// Size variable for CONMIN arrays. See CONMIN manual.
  ///** N3 = Maximum possible number of active constraints.*/
  //int N3 = xsize + 1; // + NCON;;
  ///// Size variable for CONMIN arrays. See CONMIN manual.
  ///** N4 = Maximum(N3,number of variables) */
  //int N4 = N3;
  ///// Size variable for CONMIN arrays. See CONMIN manual.
  ///** N5 = 2*(N4) */
  //int N5 = 2 * N4;
  /// Array size parameter needed in interface to CONMIN.
  //int conminSingleArray = 1;
  int NDV = static_cast<int>(ndv);
  /// CONMIN variable: Number of constraints.
  int NCON = 0;
  int N1 = NDV+2;
  int N2 = 2*NDV+NCON;
  int N3 = N2; //NDV + 2;
  int N4 = N3;
  int N5 = 2*N4;
  //printf("%d %d %d %d %d %d\n",data.xSize(),N1,N2,N3,N4,N5);
  //int N1,N2,N3,N4,N5;
  //N1 = 6;:
  //N2 = N3 = N4 = 11;
  //N5 = 22;
  /// CONMIN variable: Finite difference flag.
  int NFDG = 0;
  /// CONMIN variable: Flag to control amount of output data.
  int IPRINT = 2;
  /// CONMIN variable: Flag to specify the maximum number of iterations.
  int ITMAX = max_iter;
  /// CONMIN variable: Relative finite difference step size.
  double FDCH = 1.0e-2;
  /// CONMIN variable: Absolute finite difference step size.
  double FDCHM = 1.0e-2;
  /// CONMIN variable: Constraint thickness parameter
  /** The value of CT decreases in magnitude during optimization.*/
  double CT = -0.1;
  /// CONMIN variable: Minimum absolute value of CT used during optimization.
  double CTMIN = 0.004;
  /// CONMIN variable: Constraint thickness parameter for linear and
  /// side constraints.
  double CTL = -0.01;
  /// CONMIN variable: Minimum value of CTL used during optimization.
  double CTLMIN = 0.001;
  /// CONMIN variable: Relative convergence criterion threshold.
  /*** Threshold for the minimum relative change in the objective function. */
  double DELFUN = .001;
  /// CONMIN variable: Absolute convergence criterion threshold.
  /*** Threshold for the minimum relative change in the objective function. */
  double DABFUN = 1.0e-10;
  /// CONMIN variable: status flag for optimization
  int conminInfo = 0;
  /// Internal CONMIN array.
  /** Move direction in N-dimensional space.*/
  VecDbl S(N1,0.0);
  /// Internal CONMIN array.
  /** Temporary storage of constraint values.*/
  VecDbl G1(N2,0.0);
  /// Internal CONMIN array.
  /** Temporary storage of constraint values.*/
  VecDbl G2(N2,0.0);
  /// Internal CONMIN array.
  /** Temporary storage of constraint values.*/
  VecDbl B(N3 * N3,0.0);
  /// Internal CONMIN array.
  /** Temporary storage for use with arrays B and S.*/
  VecDbl C(N4,0.0);
  /// Internal CONMIN array.
  /** Temporary storage for use with arrays B and S.*/
  vector<int> MS1(N5,0);
  /// Internal CONMIN array.
  /** Vector of scaling parameters for design parameter values.*/
  VecDbl SCAL(N1,0.0);
  /// Internal CONMIN array.
  /** Temporary storage for analytic gradient data.*/
  VecDbl DF(N1,0.0);
  /// Internal CONMIN array.
  /** Temporary 2-D array for storage of constraint gradients.*/
  VecDbl A(N1 * N3,0.0);
  /// Internal CONMIN array.
  /** Array of flags to identify linear constraints. (not used in this
      implementation of CONMIN) */
  vector<int> ISC(N2,0);
  /// Internal CONMIN array.
  /** Array of flags to identify active and violated constraints */
  vector<int> IC(N3,0);
                                                                                  
  /// Internal CONMIN variable: 1-D search parameter.
  double ALPHAX = 0.1;
  /// Internal CONMIN variable: 1-D search parameter.
  double ABOBJ1 = 0.1;
  /// Internal CONMIN variable: mean value of push-off factor.
  double THETA = 1.0;
  /// Internal CONMIN variable: "participation coefficient".
  double PHI = 5.0;
  /// Internal CONMIN variable: side constraints parameter.
  int  NSIDE = 0;
  /// Internal CONMIN variable: scaling control parameter.
  int  NSCAL = 0;
  /// Internal CONMIN variable: estimate of 1+(max # of active constraints).
  //int  NACMX1 = N3;
  /// Internal CONMIN variable: linear objective function identifier (unused).
  int  LINOBJ = 0;
  /// Internal CONMIN variable: diminishing return criterion iteration number.
  int  ITRM = 3;
  /// Internal CONMIN variable: conjugate direction restart parameter.
  int  ICNDIR = NDV + 1;
  /// Internal CONMIN variable: internal optimization termination flag.
  int  IGOTO = 1;
  /// Internal CONMIN variable: number of active and violated constraints.
  int  NAC;
  /// Internal CONMIN variable: gradient information flag.
  int  INFOG = 0;
  /// Internal CONMIN variable: iteration count.
  int  ITER = 0;
  VecDbl query_pt(N1,1.0);
  copy(x.begin(),x.end(),query_pt.begin()); 
  
  vector<double> lower_bounds(N1,0.0);
  vector<double> upper_bounds(N1,0.0);
  lowerBounds = VecDbl(NDV,1e-3);
  upperBounds = VecDbl(NDV,1e+6);
  NSIDE = 1;
  copy(lowerBounds.begin(),lowerBounds.end(),lower_bounds.begin());
  copy(upperBounds.begin(),upperBounds.end(),upper_bounds.begin());
  cout << "Lower bounds: ";
  copy(lower_bounds.begin(),lower_bounds.end(),std::ostream_iterator<double>(cout,"\n"));
  cout << "Upper bounds: ";
  copy(upper_bounds.begin(),upper_bounds.end(),std::ostream_iterator<double>(cout,"\n"));
  VecDbl cv(N2,0.0);
  double OBJ = 0.0;
  //KrigingCPPSurface kcs(sd);
  //DABFUN = .001 * objective(x);
  do {
    //kcs.usePreComputedCorrelationVector(candidateCorrelations);
    //kcs.createModel();
    OBJ = objective(x);
    //gradient(x,DF);
    //constraints(x,A,IC,NAC,cv,CT);
    //cv[0] = con1(x);
    //cv[1] = con2(x);
    //cv[2] = con3(x);
    //cout << "OBJ: "  << OBJ << endl;
    //copy(x.begin(),x.end(),std::ostream_iterator<double>(cout,"\n"));
    //printf("OBJ: %f Constraints: %f %f %f\n",OBJ,cv[0],cv[1],cv[2]);
    //cout << "Conmin info: " << conminInfo << endl;
    //break;
    CONMIN_F77(
              &query_pt[0],&lower_bounds[0],&upper_bounds[0],
              &cv[0],
              &SCAL[0],&DF[0],&A[0],&S[0],&G1[0],&G2[0],&B[0],&C[0],
              &ISC[0],&IC[0],&MS1[0],N1,N2,N3,N4,N5,
              DELFUN,DABFUN,FDCH,FDCHM,CT,CTMIN,
              CTL,CTLMIN,ALPHAX,ABOBJ1,THETA,
              OBJ,
              NDV,NCON,NSIDE,IPRINT,NFDG,NSCAL,LINOBJ,ITMAX,
              ITRM,ICNDIR,IGOTO,NAC,conminInfo,INFOG,ITER);
    for(unsigned i = 0; i < x.size(); i++) x[i] = query_pt[i];
  printf("x.size(): %d query_pt.size(): %d N1: %d\n",x.size(),query_pt.size(),N1);
  } while (IGOTO != 0);
  final_val = OBJ;
}

double ConminKriging::objective(const VecDbl& x)
{
  const VecDbl& correlationVector = x; // this is what we're trying to optimize
  // with respect to
  MtxDbl R = KrigingModel::corrMtx(correlationVector,data); // correlation matrix
  
  VecDbl y = data.getResponses(); // vector of observed values
  // Invert the correlation matrix
  vector<int> ipvt;
  surfpack::LUFact(R,ipvt);
  //cout << R.asString() << endl;
  // Calculate the determinant before doing the inversion
  double det_corr_mtx = 1.0;
  for (unsigned i = 0; i < data.size(); i++) {
    det_corr_mtx *= R(i,i);
  }
  surfpack::inverseAfterLUFact(R,ipvt);
  VecDbl Rinv_times_y;
  surfpack::matrixVectorMult(Rinv_times_y,R,y);
  double denominator_sum = 0;
  for (unsigned i = 0; i < data.size(); i++) {
    for (unsigned j = 0; j < data.size(); j++) {
      denominator_sum += R(i,j); 
    }
  }
  double numerator_sum = surfpack::sum_vector(Rinv_times_y);
  double beta_hat = numerator_sum / denominator_sum;
  //cout << "beta_hat: " << beta_hat << endl;
  //cout << "numerator: " << numerator_sum << endl;
  //cout << "denominator: " << denominator_sum << endl;
  surfpack::vectorShift(y,beta_hat);
  surfpack::matrixVectorMult(this->rhs,R,y);
  double estVariance = surfpack::dot_product(y,this->rhs);
  //cout << "EstVariance: " << estVariance << endl;
  //cout << "det_corr_mtx: " << det_corr_mtx << endl;
  double likelihood = 
    (data.size()*log(estVariance)+log(fabs(det_corr_mtx)));
  cout << "Objective value: " << likelihood << endl;
  return likelihood;
}

VecDbl ConminKriging::gradient(const VecDbl& x)
{
  assert(false); // this function should not be called, it needs to be implemented properly first.  Frankly, I'm not sure that there is an analytical form for the derivative.
  return VecDbl(ndv,0.0);
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
  : SurfpackModelFactory(), 
  correlations(0), max_iter(50), conmin_seed(0)
{

}

KrigingModelFactory::KrigingModelFactory(const ParamMap& args)
  : SurfpackModelFactory(args), 
  correlations(0), max_iter(50), conmin_seed(0)
{

}

void KrigingModelFactory::config()
{
  SurfpackModelFactory::config();
  string strarg;
  strarg = params["correlations"];
  if (strarg != "") correlations = surfpack::toVec<double>(strarg); 
  strarg = params["max_iter"];
  if (strarg != "") {
    max_iter = std::atoi(strarg.c_str()); 
    assert(max_iter >= 0);
  }
  strarg = params["conmin_seed"];
  if (strarg != "") conmin_seed = surfpack::toVec<double>(strarg); 
  
}

SurfpackModel* KrigingModelFactory::Create(const std::string& model_string)
{
  ///\todo Be able to parse an RBF model from a string
  assert(false);
  return 0;
}


typedef std::pair<double,VecDbl> KMPair;
SurfpackModel* KrigingModelFactory::Create(const SurfData& sd)
{
  this->add("ndims",surfpack::toString(sd.xSize()));
  this->config();

  if (!correlations.empty()) {
    // don't do anything to them
  } else if (!conmin_seed.empty()) {
    correlations = conminCorrelations(sd);
  } else {
    correlations = sampleCorrelations(sd);
  }
  assert(correlations.size() == ndims);
  return new KrigingModel(sd,correlations);
}


VecDbl KrigingModelFactory::sampleCorrelations(const SurfData& sd)
{
  vector<AxesBounds::Axis> axes;
  double max_correlation = 15.0;
  for (unsigned i = 0; i < sd.xSize(); i++) {
    axes.push_back(AxesBounds::Axis(0,max_correlation));
  }
  AxesBounds ax(axes);
  KMPair best(std::numeric_limits<double>::max(),VecDbl(sd.xSize(),1.0));
  SurfData* corrs = ax.sampleMonteCarlo(5);
  
  StandardFitness sf;
  for (unsigned i = 0; i < corrs->size(); i++) {
    VecDbl candidates = (*corrs)(i);
  //cout << "Our Correlations: ";
  //copy(correlations.begin(),correlations.end(),std::ostream_iterator<double>(cout," "));
  
    KrigingModel km(sd,candidates);
    double fitness = sf(km,sd);
    cout << "Our Fitness: " << fitness << "*********************************" << endl;
    double lk = km.likelihood;
    if (lk == lk && lk > (-std::numeric_limits<double>::max()) && lk < best.first)     {
      best = KMPair(lk,candidates);
      //cout << "Update" << endl;
    }
    //} else if (lk != lk) {
    //  // NaN
    //  //cout << "Not equal to self: " << lk << endl;
    //} else if (lk >= best.first) {
    //  //cout << "Not good enough: " << lk << " " << best.first <<  endl;
    //} else if (!(lk > std::numeric_limits<double>::min())) {
    //  // -inf
    //  //cout << "Not greater than min: " << lk << endl;
    //}
  }
  copy(best.second.begin(),best.second.end(),std::ostream_iterator<double>(cout," ")); 
  //cout << " lik: " << best.first << endl;
  delete corrs;
  return best.second; 
}

VecDbl KrigingModelFactory::conminCorrelations(const SurfData& sd)
{
  ConminKriging ck(sd);
  if (conmin_seed.empty()) conmin_seed = VecDbl(sd.xSize(),1.0);
  VecDbl results = conmin_seed;
  double final;
  int max_iter = 40;
  ck.optimize(results, final, max_iter);
  return results;
}
