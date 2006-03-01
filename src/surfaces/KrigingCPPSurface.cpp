/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */
#include "surfpack_config.h"
#include "surfpack.h"
#include "SurfData.h"
#include "KrigingCPPSurface.h"
#include "SurfpackMatrix.h"

//#define PRINT_DEBUG

using namespace std;
using namespace surfpack;

#define CONMIN_F77 F77_FUNC(conmin,CONMIN)
#ifdef __cplusplus
extern "C"  /* prevent C++ name mangling */
#endif
void CONMIN_F77(double* candidate, double* lowerb, double* upperb,
		double* constraint_values,
	     double* scal, double* df, double* a, double* s, double* g1,
	     double* g2, double* b, double* c,
	     int* isc, int* ic, int* ms1,
	     int& n1, int& n2, int& n3, int& n4, int& n5,
		double& delfun, double& dabfun, double& fdch, double& fdchm,
	     double& ct, double& ctmin, double& ctl, double& ctlmin,
	     double& alphax, double& abobj1, double& theta,
	     double& obj,
	     int& numdv, int& numcon, int& nside, int& iprint, int& nfdg,
	     int& nscal, int& linobj, int& itmax, int& itrm, int& incdir,
	     int& igoto, int& nac, int& info, int& infog, int& iter);

const string KrigingCPPSurface::name = "Kriging";
//_____________________________________________________________________________
// Creation, Destruction, Initialization
//_____________________________________________________________________________

KrigingCPPSurface::KrigingCPPSurface(SurfData* sd)
  : Surface(sd), runConminFlag(true)
{
}

KrigingCPPSurface::KrigingCPPSurface(const string filename) : Surface(0)
{
  read(filename);
}

KrigingCPPSurface::~KrigingCPPSurface() 
{
}

//_____________________________________________________________________________
// Queries 
//_____________________________________________________________________________

const string KrigingCPPSurface::surfaceName() const
{
  return name;
}

unsigned KrigingCPPSurface::minPointsRequired(unsigned hypothetical_xsize) 
{ 
  return hypothetical_xsize * hypothetical_xsize;
}

unsigned KrigingCPPSurface::minPointsRequired() const
{ 
  if (this->xsize <= 0) {
    throw string(
      "Dimensionality of data needed to determine number of required samples."
    );
  } else {
    return minPointsRequired(this->xsize);
  }
}

double KrigingCPPSurface::evaluate(const std::vector<double>& x)
{ 
  vector<double> rx(sd->size()); // vector of observed values
  for (unsigned i = 0; i < sd->size(); i++) {
    rx[i] = correlation_function(correlationVector,x,(*sd)[i].X());
  }
  double result = surfpack::dot_product(rx,this->rhs);
  return this->betaHat + result;
  
}
    
//_____________________________________________________________________________
// Commands 
//_____________________________________________________________________________

void KrigingCPPSurface::setConminThetaVars(const std::vector<double>& vals)
{
  if (sd && sd->xSize() != vals.size()) {
    throw string("Dimension mismatch: conmin seed and data dimensionality");
  } else { 
    conminCorrelations = vals;
  }
}

void KrigingCPPSurface::useUniformCorrelationValue(double correlation)
{
  if (xsize == 0) {
    throw string("Must know data arity to use uniform correlation value.");
  }
  std::vector<double> vals(xsize);
  for (unsigned i = 0; i < xsize; i++) correlationVector[i] = correlation;
}
void KrigingCPPSurface::
  usePreComputedCorrelationVector(const std::vector<double>& vals)
{
  if (sd && sd->xSize() > vals.size()) {
    throw string("Dimension mismatch: correlations and data dimensionality");
  }
  correlationVector.assign(vals.begin(),vals.begin()+sd->xSize());
  runConminFlag = false;
}

void KrigingCPPSurface::build(SurfData& data)
{
  if (runConminFlag) {
    correlationVector = useConminToFindCorrelationParams();
  }
  SurfpackMatrix<double> R(data.size(),data.size(),true); // correlation matrix
  vector<double> y(data.size()); // vector of observed values
  for (unsigned i = 0; i < data.size(); i++) {
    y[i] = data.getResponse(i);
    for (unsigned j = 0; j < data.size(); j++) {
      if (i == j) {
        R[i][j] = 1.0;
      } else {
        R[i][j] = 
   	  correlation_function(correlationVector,data[i].X(),data[j].X());
      }
    }
  }
  // Invert the correlation matrix
  vector<int> ipvt;
  surfpack::LUFact(R,ipvt);
  cout << R.asString() << endl;
  // Calculate the determinant before doing the inversion
  determinantCorrMatrix = 1.0;
  for (unsigned i = 0; i < data.size(); i++) {
    determinantCorrMatrix *= R[i][i];
  }
  surfpack::inverseAfterLUFact(R,ipvt);
  vector<double> Rinv_times_y;
  surfpack::matrixVectorMult(Rinv_times_y,R,y);
  double denominator_sum = 0;
  for (unsigned i = 0; i < data.size(); i++) {
    for (unsigned j = 0; j < data.size(); j++) {
      denominator_sum += R[i][j]; 
    }
  }
  double numerator_sum = surfpack::sum_vector(Rinv_times_y);
  this->betaHat = numerator_sum / denominator_sum;
  //cout << "betaHat: " << betaHat << endl;
  //cout << "numerator: " << numerator_sum << endl;
  //cout << "denominator: " << denominator_sum << endl;
  surfpack::vectorShift(y,betaHat);
  matrixVectorMult(this->rhs,R,y);
  double estVariance = surfpack::dot_product(y,rhs);
  //cout << "EstVariance: " << estVariance << endl;
  //cout << "determinantCorrMatrix: " << determinantCorrMatrix << endl;
  this->likelihood = 
    -0.5*(sd->size()*log(estVariance)+log(abs(determinantCorrMatrix)));
}

void KrigingCPPSurface::config(const Arg& arg)
{
  vector<double> temp;
  string argname = arg.name;
  if (argname == "conmin_seed") {
    setConminThetaVars(RvalTuple::asVectorDouble(temp,arg.getRVal()->getTuple())); 
  } else if (argname == "correlations") {
    usePreComputedCorrelationVector(RvalTuple::asVectorDouble(temp,arg.getRVal()->getTuple()));
  } else if (argname == "uniform_correlation") {
    useUniformCorrelationValue(arg.getRVal()->getReal());
  } else {
    Surface::config(arg);
  }
}
/// Create a surface of the same type as 'this.'  This objects data should
/// be replaced with the dataItr passed in, but all other attributes should
/// be the same (e.g., a second-order polynomial should return another 
/// second-order polynomial.  Surfaces returned by this method can be used
/// to compute the PRESS statistic.
KrigingCPPSurface* KrigingCPPSurface::makeSimilarWithNewData(SurfData* surfData)
{
  return new KrigingCPPSurface(surfData);
}

//_____________________________________________________________________________
// Helper Functions 
//_____________________________________________________________________________

double KrigingCPPSurface::
  correlation_function(const std::vector<double>& correlations,
  const std::vector<double>& pt1, const std::vector<double>& pt2)
{
  double sum = 0.0;
  double mult_term = 0.0;
  for (unsigned i = 0; i < correlations.size(); i++) {
    mult_term = pt1[i]-pt2[i];
    sum += correlations[i]*mult_term*mult_term;
  }
  return exp(-sum);
}

std::vector<double> KrigingCPPSurface::useConminToFindCorrelationParams()
{
  //*******************CONMIN DATA**************************/

  /// Size variable for CONMIN arrays. See CONMIN manual.
  /** N1 = number of variables + 2 */
  int N1 = xsize + 2;
  /// Size variable for CONMIN arrays. See CONMIN manual.
  /** N2 = number of constraints + 2*(number of variables) */
  int N2 = xsize * 2; // + numcon;
  /// Size variable for CONMIN arrays. See CONMIN manual.
  /** N3 = Maximum possible number of active constraints.*/
  int N3 = xsize + 1; // + numcon;;
  /// Size variable for CONMIN arrays. See CONMIN manual.
  /** N4 = Maximum(N3,number of variables) */
  int N4 = N3;
  /// Size variable for CONMIN arrays. See CONMIN manual.
  /** N5 = 2*(N4) */
  int N5 = 2 * N4;
  /// Array size parameter needed in interface to CONMIN.
  //int conminSingleArray = 1;
  /// CONMIN variable: Number of constraints.
  int numcon = 0;
  /// CONMIN variable: Finite difference flag.
  int NFDG = 0;
  /// CONMIN variable: Flag to control amount of output data.
  int IPRINT = 0;  
  /// CONMIN variable: Flag to specify the maximum number of iterations.
  int ITMAX = 100;  
  /// CONMIN variable: Relative finite difference step size.
  double FDCH = 1.0e-5; 
  /// CONMIN variable: Absolute finite difference step size.
  double FDCHM = 1.0e-5;
  /// CONMIN variable: Constraint thickness parameter 
  /** The value of CT decreases in magnitude during optimization.*/
  double CT = -0.1;
  /// CONMIN variable: Minimum absolute value of CT used during optimization.
  double CTMIN = 0.001;
  /// CONMIN variable: Constraint thickness parameter for linear and
  /// side constraints.
  double CTL = -0.01;
  /// CONMIN variable: Minimum value of CTL used during optimization.
  double CTLMIN = 0.001;
  /// CONMIN variable: Relative convergence criterion threshold.
  /*** Threshold for the minimum relative change in the objective function. */
  double DELFUN = 1.0e-7; 
  /// CONMIN variable: Absolute convergence criterion threshold.
  /*** Threshold for the minimum relative change in the objective function. */
  double DABFUN = 1.0e-7;
  /// CONMIN variable: status flag for optimization
  int conminInfo = 0;
  /// Internal CONMIN array.
  /** Move direction in N-dimensional space.*/
  vector<double> S(N1);
  /// Internal CONMIN array.
  /** Temporary storage of constraint values.*/
  vector<double> G1(N2);
  /// Internal CONMIN array.
  /** Temporary storage of constraint values.*/
  vector<double> G2(N2);
  /// Internal CONMIN array.
  /** Temporary storage for computations involving array S.*/
  vector<double> B(N3 * N3);
  /// Internal CONMIN array.
  /** Temporary storage for use with arrays B and S.*/
  vector<double> C(N4);
  /// Internal CONMIN array.
  /** Temporary storage for use with arrays B and S.*/
  vector<int> MS1(N5);
  /// Internal CONMIN array.
  /** Vector of scaling parameters for design parameter values.*/
  vector<double> SCAL(N1);
  /// Internal CONMIN array.
  /** Temporary storage for analytic gradient data.*/
  vector<double> DF(N1);
  /// Internal CONMIN array.
  /** Temporary 2-D array for storage of constraint gradients.*/
  vector<double> A(N1 * N3);
  /// Internal CONMIN array.
  /** Array of flags to identify linear constraints. (not used in this
      implementation of CONMIN) */
  vector<int> ISC(N2); 
  /// Internal CONMIN array.
  /** Array of flags to identify active and violated constraints */
  vector<int> IC(N3);

  /// Internal CONMIN variable: 1-D search parameter.
  double ALPHAX = 0.1;
  /// Internal CONMIN variable: 1-D search parameter.
  double ABOBJ1 = 0.1;
  /// Internal CONMIN variable: mean value of push-off factor.
  double THETA = 1.0;
  /// Internal CONMIN variable: "participation coefficient".
  //double PHI = 5.0;
  /// Internal CONMIN variable: side constraints parameter.
  int  NSIDE = 1;
  /// Internal CONMIN variable: scaling control parameter.
  int  NSCAL = 0;
  /// Internal CONMIN variable: estimate of 1+(max # of active constraints).
  //int  NACMX1 = N3;
  /// Internal CONMIN variable: linear objective function identifier (unused).
  int  LINOBJ = 0;
  /// Internal CONMIN variable: diminishing return criterion iteration number.
  int  ITRM = 3;
  /// Internal CONMIN variable: conjugate direction restart parameter.
  int  ICNDIR = xsize + 1;
  /// Internal CONMIN variable: internal optimization termination flag.
  int  IGOTO = 1;
  /// Internal CONMIN variable: number of active and violated constraints.
  int  NAC;
  /// Internal CONMIN variable: gradient information flag.
  int  INFOG;
  /// Internal CONMIN variable: iteration count.
  int  ITER = 0;
  
  int numdv = static_cast<int>(xsize);
  vector<double> candidateCorrelations(N1,1.0);
  vector<double> lower_bounds(N1,0.0);
  vector<double> upper_bounds(N1,1.0e3);
  vector<double> constraintVector(N2);
  double OBJ = 0.0;
  KrigingCPPSurface kcs(sd);
  do {
    kcs.usePreComputedCorrelationVector(candidateCorrelations);
    kcs.createModel();
    OBJ = kcs.likelihood;
    cout << "OBJ: "  << OBJ << endl;
    CONMIN_F77(
              &candidateCorrelations[0],&lower_bounds[0],&upper_bounds[0],
              &constraintVector[0],
              &SCAL[0],&DF[0],&A[0],&S[0],&G1[0],&G2[0],&B[0],&C[0],
              &ISC[0],&IC[0],&MS1[0],N1,N2,N3,N4,N5,
              DELFUN,DABFUN,FDCH,FDCHM,CT,CTMIN,
              CTL,CTLMIN,ALPHAX,ABOBJ1,THETA,
              OBJ,
              numdv,numcon,NSIDE,IPRINT,NFDG,NSCAL,LINOBJ,ITMAX,
              ITRM,ICNDIR,IGOTO,NAC,conminInfo,INFOG,ITER);
  } while (IGOTO != 0);
  return candidateCorrelations;
}

double KrigingCPPSurface::likelihoodEstimation()
{
  assert(sd);
  double mse = this->goodnessOfFit(string("mean_squared"),sd);
  return -0.5*(sd->size()*log(mse)+log(abs(determinantCorrMatrix)));
  
}
//_____________________________________________________________________________
// I/O 
//_____________________________________________________________________________

void KrigingCPPSurface::writeBinary(ostream& os)
{
//  int i;
//  unsigned j;
//  os.write((char*)(&xsize),sizeof(int));
//  os.write((char*)(&numsamp),sizeof(int));
//  for(i=0;i<numsamp;i++) { 
//    os.write(reinterpret_cast<char*>(&xMatrix[i]),sizeof(xMatrix[i]));
//  }
//  for(i=0;i<numsamp;i++) { 
//    os.write(reinterpret_cast<char*>(&rXhatVector[i]),sizeof(rXhatVector[i]));
//  }
//  os.write(reinterpret_cast<char*>(&betaHat),sizeof(betaHat));
//  for(i=0;i<numsamp;i++) { 
//    os.write(reinterpret_cast<char*>(&rhsTermsVector[i]),sizeof(rhsTermsVector[i]));
//  }
//  for(j=0;j<xsize;j++) { 
//    os.write(reinterpret_cast<char*>(&thetaVector[j]),sizeof(thetaVector[j]));
//  }
}

void KrigingCPPSurface::writeText(ostream& os)
{
//  ios::fmtflags old_flags = os.flags();
//  unsigned old_precision = os.precision(surfpack::output_precision);
//  os.setf(ios::scientific);
//  int i;
//  unsigned j;
//  os << numsamp << " number of data points" << endl;
//  os << xsize << " number of input variables" << endl;
//  for(i=0;i<numsamp;i++) { 
//    os << xMatrix[i] << " xMatrix[" << i << "]" << endl; 
//  }
//  for(i=0;i<numsamp;i++) { 
//    os <<  rXhatVector[i] << " rXhatVector[" << i << "]" << endl; 
//  }
//  os << betaHat << " betaHat" << endl; 
//  for(i=0;i<numsamp;i++) { 
//    os <<  rhsTermsVector[i] << " rhsTermsVector[" << i << "]" << endl; 
//  }
//  for(j=0;j<xsize;j++) { 
//    os << thetaVector[j] << " thetaVector[" << j << "]" << endl; 
//  }
//  os.flags(old_flags);
//  os.precision(old_precision);
}

void KrigingCPPSurface::readBinary(istream& is)
{
//  int i;
//  unsigned j;
//  is.read((char*)(&xsize),sizeof(int));
//  is.read((char*)(&numsamp),sizeof(int));
//  initialize();
//  for(i=0;i<numsamp;i++) { 
//    is.read(reinterpret_cast<char*>(&xMatrix[i]),sizeof(xMatrix[i]));
//  }
//  for(i=0;i<numsamp;i++) { 
//    is.read(reinterpret_cast<char*>(&rXhatVector[i]),sizeof(rXhatVector[i]));
//  }
//  is.read(reinterpret_cast<char*>(&betaHat),sizeof(betaHat));
//  for(i=0;i<numsamp;i++) { 
//    is.read(reinterpret_cast<char*>(&rhsTermsVector[i]),sizeof(rhsTermsVector[i]));
//  }
//  for(j=0;j<xsize;j++) { 
//    is.read(reinterpret_cast<char*>(&thetaVector[j]),sizeof(thetaVector[j]));
//  }
}

void KrigingCPPSurface::readText(istream& is)
{
//  string sline;
//  istringstream streamline;
//  int i;
//  unsigned j;
//  getline(is,sline); streamline.str(sline);
//  streamline >> numsamp;
//  getline(is,sline); streamline.str(sline);
//  streamline >> xsize;
//  initialize();
//  for(i=0;i<numsamp;i++) { 
//    getline(is,sline); streamline.str(sline);
//    streamline >> xMatrix[i];
//  }
//  for(i=0;i<numsamp;i++) { 
//    getline(is,sline); streamline.str(sline);
//    streamline >> rXhatVector[i];
//  }
//  getline(is,sline); streamline.str(sline);
//  streamline >> betaHat;
//  for(i=0;i<numsamp;i++) { 
//    getline(is,sline); streamline.str(sline);
//    streamline >> rhsTermsVector[i];
//  }
//  for(i=0;i<xsize;i++) { 
//    getline(is,sline); streamline.str(sline);
//    streamline >> thetaVector[i];
//  }
//
}

//_____________________________________________________________________________
// Testing 
//_____________________________________________________________________________

#ifdef __TESTING_MODE__
#endif
