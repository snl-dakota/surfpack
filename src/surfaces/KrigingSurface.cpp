#include <cmath>
#include <cfloat>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

#include "SurfPoint.h"
#include "SurfData.h"
#include "Surface.h"
#include "SkipSurfDataIterator.h"
#include "SurfDataIterator.h"
#include "KrigingSurface.h"



#ifndef RS6K
#define callconmin callconmin_
#define krigmodel  krigmodel_
#endif

//#define PRINT_DEBUG

// Prototypes for the DAKOTA-to-CONMIN interface subroutine "callconmin"
// and for the kriging model evaluation subroutine "krigmodel".
extern "C" void callconmin(double&, double&, double&, double&, double&, double&,
			   double&, double&, double&, double&, double&, double&,
			   int&, int&, int&, int&, int&, int&, int&, int&,
			   double&, double&, double&, double&, double&, double&,
			   double&, double&, double&, double&, double&, double&,
			   int&, int&, int&, int&, int&, int&, int&, int&,
			   int&, int&, int&, int&, int&, int&, int&,
			   int&, int&, int&,
			   double&, double&, double&, double&, double&, double&, int&,
			   double&, double&, double&, double&, double&, double&, 
			   double&, double&, double&, int&, int&, int&);

extern "C" void krigmodel(int&, int&, int&,
			  int&, double&, double&, double&, const double&,
			  double&, double&, double&, double&, int&, 
			  double&, double&, double&, double&, double&, 
			  double&, double&, double&, double&, int&, int&);

using namespace std;

const string KrigingSurface::name = "Kriging";
//_____________________________________________________________________________
// Creation, Destruction, Initialization
//_____________________________________________________________________________

KrigingSurface::KrigingSurface(AbstractSurfDataIterator* dataItr) 
  : Surface(dataItr), needsCleanup(false), runConminFlag(true)
{

#ifdef __TESTING_MODE__
  constructCount++;
#endif
  build();
}
KrigingSurface::KrigingSurface(SurfData& sd, unsigned responseIndex)
  : Surface(0), needsCleanup(false), runConminFlag(true)
{
#ifdef __TESTING_MODE__
  constructCount++;
#endif
  dataItr = new SurfDataIterator(sd, responseIndex);
  numdv = numsamp = 0;
  build();
}

KrigingSurface::KrigingSurface(const string filename) : Surface(0)
{
#ifdef __TESTING_MODE__
  constructCount++;
#endif
  read(filename);
}

KrigingSurface::~KrigingSurface() 
{
#ifdef __TESTING_MODE__
  destructCount++;
#endif
   if (needsCleanup) {
      cleanup();
   }
}

void KrigingSurface::initialize()
{
  if (numdv <= 0 || numsamp <= 0) {
    cerr << "Trying to create Kriging Surface with bad or non-existent surfdata" << endl;
  }
  // CONMIN parameters and array allocation
  NFDG   = 0;       // default finite difference flag
  IPRINT = 0;       // default flag to control amount of output info
  ITMAX  = 100;     // default max number of iterations
  FDCH   =  1.0e-5; // default relative finite difference step size
  FDCHM  =  1.0e-5; // default absolute finite difference step size
  CT     = -0.1;    // default constraint "thickness tolerance 
                    // (for determining active/inactive constraint status)
  CTMIN  =  0.001;  // default absolute constraint tolerance 
                    // (note: the CONMIN manual default is 0.004)
  CTL    = -0.01;   // default side constraint thickness tolerance (see CT)
  CTLMIN =  0.001;  // default absolute side constraint tolerance
  DELFUN =  1.0e-7; // default minimum relative change in the objective
                    // function needed for convergence
  DABFUN =  1.0e-7; // default minimum absolute change in the objective
                    // function needed for convergence
  conminInfo = 0;   // must be set to 0 before calling CONMIN
  NSIDE  = 1;       // flag to identify existence of side constraints in 
                    // optimization subproblem to solve for kriging 
                    // correlation parameters

  numcon = 0;     // number of constraints
  N1     = numdv + 2;
  N2     = numdv*2 + numcon;
  N3     = 1 + numcon + numdv;
  N4     = (N3 >= numdv) ? N3 : numdv;
  N5     = 2*N4;
  S      = new double[N1];
  G1     = new double[N2];
  G2     = new double[N2];
  B      = new double[N3*N3];
  C      = new double[N4];
  MS1    = new int[N5];
  SCAL   = new double[N1];
  DF     = new double[N1];
  A      = new double[N1*N3];
  ISC    = new int[N2];
  IC     = new int[N3];
  conminThetaVars      = new double[N1];
  conminThetaLowerBnds = new double[N1];  
  conminThetaUpperBnds = new double[N1];
  
  ICNDIR   =  numdv+1;  // conjugate direction restart parameter
  NSCAL    =  0;  // (unused) scaling flag
  NACMX1   =  N3; // estimate of 1+(max # of active constraints)
  LINOBJ   =  0;  // (unused) set to 1 if obj fcn is linear
  ITRM     =  3;  // number of consecutive iterations of less than DELFUN
                  // or DABFUN progress before optimization terminated
  THETA    =  1.0;// mean value of "push-off factor" in mtd of feasible
                  // directions
  PHI      =  5.0;// "participation coefficient" to force an infeasible
                  // design to be feasible
  ALPHAX   = 0.1; // 1-D search fractional change parameter
  ABOBJ1   = 0.1; // 1-D search fractional change parameter for 1st step

  // initialize CONMIN arrays IC and ISC
  int i;
  for( i=0; i<numcon; i++ ) {
    IC[i]  = 0;
    ISC[i] = 0;
  }

  // Kriging parameters and array allocation
  iFlag             = 1;
  betaHat           = 0.0;
  maxLikelihoodEst  = 0.0;
  numNewPts         = 1;
  numSampQuad       = 4*numsamp;
  conminSingleArray = 1;
  xNewVector        = new double[numdv*numNewPts];
  thetaVector       = new double[numdv];
  xArray 	    = new double[numdv];
  xMatrix           = new double[numdv*numsamp];
  yValueVector      = new double[numsamp];
  rhsTermsVector    = new double[numsamp];
  constraintVector  = new double[conminSingleArray];
  thetaLoBndVector  = new double[numdv];
  thetaUpBndVector  = new double[numdv];
  iPivotVector      = new int[numsamp];
  correlationMatrix = new double[numsamp*numsamp];
  invcorrelMatrix   = new double[numsamp*numsamp];
  fValueVector      = new double[numsamp];
  fRinvVector       = new double[numsamp];
  yfbVector         = new double[numsamp];
  yfbRinvVector     = new double[numsamp];
  rXhatVector       = new double[numsamp];
  workVector        = new double[numsamp];
  workVectorQuad    = new double[4*numsamp];
  iworkVector       = new int[numsamp];
  yNewVector        = new double[numNewPts];

  // Initialize 'thetaVector' which is the vector of correlation 
  // parameters for the kriging model -- maximum liklihood estimation
  // is used to find the optimal values of thetaVector in the kriging code.
  // Also initialize the lower and upper bounds on theta.
  /*for (i=0; i<numdv; i++) {
    thetaVector[i]      = theta_input[i];
    thetaLoBndVector[i] = theta_low_bnd_input[i];
    thetaUpBndVector[i] = theta_upp_bnd_input[i];
  }*/

  // copy the "theta" vectors to the "conmin_theta" -- this is needed 
  // because CONMIN uses array sizes that are larger than the number
  // of design variables (array length = N1, and N1 = numdv+2)
  for (i=0;i<N1;i++) {
    conminThetaVars[i]      = 0.0;
    conminThetaLowerBnds[i] = -DBL_MAX;
    conminThetaUpperBnds[i] =  DBL_MAX;
  }
  for (i=0;i<numdv;i++) {
    conminThetaVars[i]      = 10.0 ;
    conminThetaLowerBnds[i] = 1.e-3; 
    conminThetaUpperBnds[i] = 1.e+16; 
  }

  // now that all this memory has been allocated, we need to set a flag
  // so that the object will know to delete it
  needsCleanup = true;
}

void KrigingSurface::cleanup()
{
  // clean up CONMIN array allocation
  delete [] conminThetaVars; conminThetaVars = 0;
  delete [] conminThetaLowerBnds; conminThetaLowerBnds = 0;
  delete [] conminThetaUpperBnds; conminThetaUpperBnds = 0;
  delete [] S; S = 0;
  delete [] G1; G1 = 0;
  delete [] G2; G2 = 0;
  delete [] B; B = 0;
  delete [] C; C = 0;
  delete [] MS1; MS1 = 0;
  delete [] SCAL; SCAL = 0;
  delete [] DF; DF = 0;
  delete [] A; A = 0;
  delete [] ISC; ISC = 0;
  delete [] IC; IC = 0;

  // clean up kriging array allocation
  delete [] xNewVector ; xNewVector = 0;
  delete [] thetaVector ; thetaVector = 0;
  delete [] xArray; xArray = 0;
  delete [] xMatrix ; xMatrix = 0;
  delete [] yValueVector ; yValueVector = 0;
  delete [] rhsTermsVector ; rhsTermsVector = 0;
  delete [] constraintVector ; constraintVector = 0;
  delete [] thetaLoBndVector; thetaLoBndVector = 0;
  delete [] thetaUpBndVector; thetaUpBndVector = 0;
  delete [] iPivotVector ; iPivotVector = 0;
  delete [] correlationMatrix ; correlationMatrix = 0;
  delete [] invcorrelMatrix ; invcorrelMatrix = 0;
  delete [] fValueVector ; fValueVector = 0;
  delete [] fRinvVector ; fRinvVector = 0;
  delete [] yfbVector ; yfbVector = 0;
  delete [] yfbRinvVector ; yfbRinvVector = 0;
  delete [] rXhatVector ; rXhatVector = 0;
  delete [] workVector ; workVector = 0;
  delete [] workVectorQuad ; workVectorQuad = 0;
  delete [] iworkVector ; iworkVector = 0;
  delete [] yNewVector ; yNewVector = 0;
}

//_____________________________________________________________________________
// Queries 
//_____________________________________________________________________________

const string KrigingSurface::surfaceName() const
{
  return name;
}

unsigned KrigingSurface::minPointsRequired(unsigned xsize) 
{ 
  return xsize + 1;
}

unsigned KrigingSurface::minPointsRequired() const
{ 
  return dataItr->xSize(); 
}

double KrigingSurface::evaluate(const std::vector<double>& x)
{ 
  //double* xArray = new double[dim];
  if (x.size() != static_cast<unsigned>(dim)) {
    cerr << "Wrong number of dimensions" << endl;
  } else {
    for (int i = 0; i < dim; i++) {
      xArray[i] = x[i];
    }
  }
  yNewVector[0] =  0.0;
  iFlag = 2; // value = 2 for evaluating a previously computed kriging model
#ifdef PRINT_DEBUG
  ostream& os = cout;
  os << "Before call to krigmodel in calculate" << endl;
  printKrigModelVariables(os);
#endif
  //cout << "Calculate_________________" << endl;
  //printKrigEvalVars(cout);
  krigmodel(dim,pts,numNewPts,
            iFlag,thetaVector[0],xMatrix[0],yValueVector[0],xArray[0],
            yNewVector[0],betaHat,rhsTermsVector[0],maxLikelihoodEst,
            iPivotVector[0],correlationMatrix[0],invcorrelMatrix[0],
            fValueVector[0],fRinvVector[0],yfbVector[0],yfbRinvVector[0],
            rXhatVector[0],workVector[0],workVectorQuad[0],iworkVector[0],
            numSampQuad);
#ifdef PRINT_DEBUG
  os << "After call to krigmodel in calculate" << endl;
  printKrigModelVariables(os);
#endif
  return yNewVector[0];
}
    
//_____________________________________________________________________________
// Commands 
//_____________________________________________________________________________

void KrigingSurface::usePreComputedCorrelationVector(std::vector<double> vals)
{
  for (unsigned i = 0; i < vals.size(); i++) {
	  thetaVector[i] = vals[i];
  }
}

void KrigingSurface::build()
{
  double* saveTheta;

  // code from Surface::build that needs to be executed before
  // initialize
  // this could be cleaner

  if (!acceptableData())
  {
    cerr << "Unacceptable data.  Cannot build Kriging surface." << endl;
  } else {
    numdv = dataItr->xSize();
    numsamp = dataItr->elementCount();
    
if (!runConminFlag) {
      saveTheta = new double[numdv];
      memcpy(saveTheta,thetaVector,sizeof(double)*numdv);
    }

    if (needsCleanup) {
      cleanup();
    }
    initialize();

    if (!runConminFlag) {
      memcpy(thetaVector,saveTheta,sizeof(double)*numdv);
      delete saveTheta;
      saveTheta = 0;
    }
    buildModel();
    valid = true;
    originalData = true;
  }
}

/// Create a surface of the same type as 'this.'  This objects data should
/// be replaced with the dataItr passed in, but all other attributes should
/// be the same (e.g., a second-order polynomial should return another 
/// second-order polynomial.  Surfaces returned by this method can be used
/// to compute the PRESS statistic.
KrigingSurface* KrigingSurface::makeSimilarWithNewData
  (AbstractSurfDataIterator* dataItr)
{
  return new KrigingSurface(dataItr);
}
//_____________________________________________________________________________
// I/O 
//_____________________________________________________________________________

void KrigingSurface::writeBinary(ostream& os)
{
  int i;
  //int tempID = Surface::krigingSurfaceID;
  
  unsigned nameSize = name.size();
  os.write(reinterpret_cast<char*>(&nameSize),sizeof(nameSize));
  os.write(name.c_str(),nameSize);
  os.write((char*)(&numdv),sizeof(int));
  os.write((char*)(&numsamp),sizeof(int));
  for(i=0;i<numdv;i++) { os.write((char*)&thetaVector[i],sizeof(double)); }
  for(i=0;i<N1;i++) { os.write((char*)&conminThetaVars[i],sizeof(double)); }
  for(i=0;i<N1;i++) { os.write((char*)&conminThetaLowerBnds[i],sizeof(double)); }
  for(i=0;i<N1;i++) { os.write((char*)&conminThetaUpperBnds[i],sizeof(double)); }
  for(i=0;i<conminSingleArray;i++) { os.write((char*)&constraintVector[i],sizeof(double)); }
  for(i=0;i<N1;i++) { os.write((char*)&SCAL[i],sizeof(double)); }
  for(i=0;i<N1;i++) { os.write((char*)&DF[i],sizeof(double)); }
  for(i=0;i<N1*N3;i++) { os.write((char*)&A[i],sizeof(double)); }
  for(i=0;i<N1;i++) { os.write((char*)&S[i],sizeof(double)); }
  for(i=0;i<N2;i++) { os.write((char*)&G1[i],sizeof(double)); }
  for(i=0;i<N2;i++) { os.write((char*)&G2[i],sizeof(double)); }
  for(i=0;i<N3*N3;i++) { os.write((char*)&B[i],sizeof(double)); }
  for(i=0;i<N4;i++) { os.write((char*)&C[i],sizeof(double)); }
  for(i=0;i<N2;i++) { os.write((char*)&ISC[i],sizeof(int)); }
  for(i=0;i<N3;i++) { os.write((char*)&IC[i],sizeof(int)); }
  for(i=0;i<N5;i++) { os.write((char*)&MS1[i],sizeof(int)); }
  os.write((char*)&N1,sizeof(int));
  os.write((char*)&N2,sizeof(int));
  os.write((char*)&N3,sizeof(int));
  os.write((char*)&N4,sizeof(int));
  os.write((char*)&N5,sizeof(int));
  os.write((char*)&N5,sizeof(int));
  os.write((char*)&DELFUN,sizeof(double));
  os.write((char*)&DABFUN,sizeof(double));
  os.write((char*)&FDCH,sizeof(double));
  os.write((char*)&FDCHM,sizeof(double));
  os.write((char*)&CT,sizeof(double));
  os.write((char*)&CTMIN,sizeof(double));
  os.write((char*)&CTL,sizeof(double));
  os.write((char*)&CTLMIN,sizeof(double));
  os.write((char*)&ALPHAX,sizeof(double));
  os.write((char*)&ABOBJ1,sizeof(double));
  os.write((char*)&THETA,sizeof(double));
  os.write((char*)&maxLikelihoodEst,sizeof(double));
  os.write((char*)&dim,sizeof(int));
  os.write((char*)&numcon,sizeof(int));
  os.write((char*)&NSIDE,sizeof(int));
  os.write((char*)&IPRINT,sizeof(int));
  os.write((char*)&NFDG,sizeof(int));
  os.write((char*)&NSCAL,sizeof(int));
  os.write((char*)&LINOBJ,sizeof(int));
  os.write((char*)&ITMAX,sizeof(int));
  os.write((char*)&ITRM,sizeof(int));
  os.write((char*)&ICNDIR,sizeof(int));
  os.write((char*)&IGOTO,sizeof(int));
  os.write((char*)&NAC,sizeof(int));
  os.write((char*)&conminInfo,sizeof(int));
  os.write((char*)&INFOG,sizeof(int));
  os.write((char*)&ITER,sizeof(int));
  os.write((char*)&pts,sizeof(int));
  os.write((char*)&numNewPts,sizeof(int));
  os.write((char*)&iFlag,sizeof(int));
  for(i=0;i<dim*pts;i++) { os.write((char*)&xMatrix[i],sizeof(double)); }
  for(i=0;i<pts;i++) { os.write((char*)& yValueVector[i],sizeof(double)); }
  for(i=0;i<dim*numNewPts;i++) { os.write((char*)& xNewVector[i] ,sizeof(double)); }
  for(i=0;i<numNewPts;i++) { os.write((char*)& yNewVector[i],sizeof(double)); }
  os.write((char*)&betaHat,sizeof(double)); 
  for(i=0;i<pts;i++) { os.write((char*)& rhsTermsVector[i],sizeof(double)); }
  for(i=0;i<pts;i++) { os.write((char*)& iPivotVector[i],sizeof(int)); }
  for(i=0;i<pts*pts;i++) { os.write((char*)& correlationMatrix[i],sizeof(double)); }
  for(i=0;i<pts*pts;i++) { os.write((char*)& invcorrelMatrix[i],sizeof(double)); }
  for(i=0;i<pts;i++) { os.write((char*)& fValueVector[i],sizeof(double)); }
  for(i=0;i<pts;i++) { os.write((char*)& fRinvVector[i],sizeof(double)); }
  for(i=0;i<pts;i++) { os.write((char*)& yfbVector[i],sizeof(double)); }
  for(i=0;i<pts;i++) { os.write((char*)& yfbRinvVector[i],sizeof(double)); }
  for(i=0;i<pts;i++) { os.write((char*)& rXhatVector[i],sizeof(double)); }
  for(i=0;i<pts;i++) { os.write((char*)& workVector[i],sizeof(double)); }
  for(i=0;i<4*pts;i++) { os.write((char*)& workVectorQuad[i],sizeof(double)); }
  for(i=0;i<pts;i++) { os.write((char*)& iworkVector[i],sizeof(int)); }
  //surfData->write(os, true); // binary write
}

void KrigingSurface::writeText(ostream& os)
{

}

void KrigingSurface::readBinary(istream& is)
{
  unsigned nameSize;
  is.read(reinterpret_cast<char*>(&nameSize),sizeof(nameSize));
  char* surfaceType = new char[nameSize+1];
  is.read(surfaceType,nameSize);
  surfaceType[nameSize] = '\0';
  string nameInFile(surfaceType);
  delete [] surfaceType;
  if (nameInFile != name) {
    cerr << "Surface name in file is not 'Kriging'." << endl;
    cerr << "Cannot build surface." << endl;
    return;
  }
  int i;
  //int tempID;
  is.read((char*)(&numdv),sizeof(int));
  dim = numdv;
  is.read((char*)(&numsamp),sizeof(int));
  pts = numsamp;
  initialize();
  for(i=0;i<numdv;i++) { is.read((char*)&thetaVector[i],sizeof(double)); }
  for(i=0;i<N1;i++) { is.read((char*)&conminThetaVars[i],sizeof(double)); }
  for(i=0;i<N1;i++) { is.read((char*)&conminThetaLowerBnds[i],sizeof(double)); }
  for(i=0;i<N1;i++) { is.read((char*)&conminThetaUpperBnds[i],sizeof(double)); }
  for(i=0;i<conminSingleArray;i++) { is.read((char*)&constraintVector[i],sizeof(double)); }
  for(i=0;i<N1;i++) { is.read((char*)&SCAL[i],sizeof(double)); }
  for(i=0;i<N1;i++) { is.read((char*)&DF[i],sizeof(double)); }
  for(i=0;i<N1*N3;i++) { is.read((char*)&A[i],sizeof(double)); }
  for(i=0;i<N1;i++) { is.read((char*)&S[i],sizeof(double)); }
  for(i=0;i<N2;i++) { is.read((char*)&G1[i],sizeof(double)); }
  for(i=0;i<N2;i++) { is.read((char*)&G2[i],sizeof(double)); }
  for(i=0;i<N3*N3;i++) { is.read((char*)&B[i],sizeof(double)); }
  for(i=0;i<N4;i++) { is.read((char*)&C[i],sizeof(double)); }
  for(i=0;i<N2;i++) { is.read((char*)&ISC[i],sizeof(int)); }
  for(i=0;i<N3;i++) { is.read((char*)&IC[i],sizeof(int)); }
  for(i=0;i<N5;i++) { is.read((char*)&MS1[i],sizeof(int)); }
  is.read((char*)&N1,sizeof(int));
  is.read((char*)&N2,sizeof(int));
  is.read((char*)&N3,sizeof(int));
  is.read((char*)&N4,sizeof(int));
  is.read((char*)&N5,sizeof(int));
  is.read((char*)&N5,sizeof(int));
  is.read((char*)&DELFUN,sizeof(double));
  is.read((char*)&DABFUN,sizeof(double));
  is.read((char*)&FDCH,sizeof(double));
  is.read((char*)&FDCHM,sizeof(double));
  is.read((char*)&CT,sizeof(double));
  is.read((char*)&CTMIN,sizeof(double));
  is.read((char*)&CTL,sizeof(double));
  is.read((char*)&CTLMIN,sizeof(double));
  is.read((char*)&ALPHAX,sizeof(double));
  is.read((char*)&ABOBJ1,sizeof(double));
  is.read((char*)&THETA,sizeof(double));
  is.read((char*)&maxLikelihoodEst,sizeof(double));
  is.read((char*)&dim,sizeof(int));
  is.read((char*)&numcon,sizeof(int));
  is.read((char*)&NSIDE,sizeof(int));
  is.read((char*)&IPRINT,sizeof(int));
  is.read((char*)&NFDG,sizeof(int));
  is.read((char*)&NSCAL,sizeof(int));
  is.read((char*)&LINOBJ,sizeof(int));
  is.read((char*)&ITMAX,sizeof(int));
  is.read((char*)&ITRM,sizeof(int));
  is.read((char*)&ICNDIR,sizeof(int));
  is.read((char*)&IGOTO,sizeof(int));
  is.read((char*)&NAC,sizeof(int));
  is.read((char*)&conminInfo,sizeof(int));
  is.read((char*)&INFOG,sizeof(int));
  is.read((char*)&ITER,sizeof(int));
  is.read((char*)&pts,sizeof(int));
  is.read((char*)&numNewPts,sizeof(int));
  is.read((char*)&iFlag,sizeof(int));
  for(i=0;i<dim*pts;i++) { is.read((char*)&xMatrix[i],sizeof(double)); }
  for(i=0;i<pts;i++) { is.read((char*)& yValueVector[i],sizeof(double)); }
  for(i=0;i<dim*numNewPts;i++) { is.read((char*)& xNewVector[i] ,sizeof(double)); }
  for(i=0;i<numNewPts;i++) { is.read((char*)& yNewVector[i],sizeof(double)); }
  is.read((char*)&betaHat,sizeof(double)); 
  for(i=0;i<pts;i++) { is.read((char*)& rhsTermsVector[i],sizeof(double)); }
  for(i=0;i<pts;i++) { is.read((char*)& iPivotVector[i],sizeof(int)); }
  for(i=0;i<pts*pts;i++) { is.read((char*)& correlationMatrix[i],sizeof(double)); }
  for(i=0;i<pts*pts;i++) { is.read((char*)& invcorrelMatrix[i],sizeof(double)); }
  for(i=0;i<pts;i++) { is.read((char*)& fValueVector[i],sizeof(double)); }
  for(i=0;i<pts;i++) { is.read((char*)& fRinvVector[i],sizeof(double)); }
  for(i=0;i<pts;i++) { is.read((char*)& yfbVector[i],sizeof(double)); }
  for(i=0;i<pts;i++) { is.read((char*)& yfbRinvVector[i],sizeof(double)); }
  for(i=0;i<pts;i++) { is.read((char*)& rXhatVector[i],sizeof(double)); }
  for(i=0;i<pts;i++) { is.read((char*)& workVector[i],sizeof(double)); }
  for(i=0;i<4*pts;i++) { is.read((char*)& workVectorQuad[i],sizeof(double)); }
  for(i=0;i<pts;i++) { is.read((char*)& iworkVector[i],sizeof(int)); }
  //delete surfData;
  //surfData = new SurfData;
  //surfData->read(is,true); // binary read
}

void KrigingSurface::readText(istream& is)
{

}

void KrigingSurface::printKrigModelVariables(ostream& os)
{
  int i;
  os << "After call to krigmodel in modelbuild" << endl;
  os << "dim: " << dim << endl;
  os << "pts: " << pts  << endl;
  os << "numNewPts: " << numNewPts << endl;
  os << "iFlag" << iFlag << endl;
  os << "thetaVector" << endl; for(i=0;i<dim;i++) { os << thetaVector[i] << endl; }
  os << "xMatrix" << endl; for(i=0;i<pts;i++) { os << xMatrix[i] << endl; }
  os << "yValueVector" << endl; for(i=0;i<pts;i++) { os <<  yValueVector[i] << endl; }
  os << "xNewVector" << endl; for(i=0;i<pts;i++) { os <<  xNewVector[i]  << endl; }
  os << "yNewVector" << endl; for(i=0;i<pts;i++) { os <<  yNewVector[i] << endl; }
  os << "betaHat: " <<  betaHat << endl; 
  os << "rhsTermsVector" << endl; for(i=0;i<pts;i++) { os <<  rhsTermsVector[i] << endl; }
  os << "maxLikelihoodEst: " << maxLikelihoodEst << endl;
  os << "iPivotVector" << endl; for(i=0;i<pts;i++) { os <<  iPivotVector[i] << endl; }
  os << "correlationMatrix" << endl; for(i=0;i<pts*pts;i++) { os <<  correlationMatrix[i] << endl; }
  os << "invcorrelMatrix" << endl; for(i=0;i<pts*pts;i++) { os <<  invcorrelMatrix[i] << endl; }
  os << "fValueVector" << endl; for(i=0;i<pts;i++) { os <<  fValueVector[i] << endl; }
  os << "fRinvVector" << endl; for(i=0;i<pts;i++) { os <<  fRinvVector[i] << endl; }
  os << "yfbVector" << endl; for(i=0;i<pts;i++) { os <<  yfbVector[i] << endl; }
  os << "yfbRinvVector" << endl; for(i=0;i<pts;i++) { os <<  yfbRinvVector[i] << endl; }
  os << "rXhatVector" << endl; for(i=0;i<pts;i++) { os <<  rXhatVector[i] << endl; }
  os << "workVector" << endl; for(i=0;i<pts;i++) { os <<  workVector[i] << endl; }
  os << "workVectorQuad" << endl; for(i=0;i<4*pts;i++) { os <<  workVectorQuad[i] << endl; }
  os << "iworkVector" << endl; for(i=0;i<pts;i++) { os <<  iworkVector[i] << endl; }
  os << "numSampQuad: " << pts  << endl;
}

void KrigingSurface::printKrigEvalVars(ostream& os)
{
  int i;
  os << "xMatrix" << endl; for(i=0;i<pts;i++) { os << xMatrix[i] << endl; }
  os << "rXhatVector" << endl; for(i=0;i<pts;i++) { os <<  rXhatVector[i] << endl; }
  os << "betaHat: " <<  betaHat << endl; 
  os << "rhsTermsVector" << endl; for(i=0;i<pts;i++) { os <<  rhsTermsVector[i] << endl; }
  os << "pts: " << pts  << endl;
  os << "dim: " << dim << endl;
  os << "numNewPts: " << numNewPts << endl;
  os << "thetaVector: " << endl; for(i=0;i<numdv;i++) { os << thetaVector[i] << endl; }
  os << "yNewVector" << endl; for(i=0;i<pts;i++) { os <<  yNewVector[i] << endl; }
}

void KrigingSurface::printConminVariables(ostream& os)
{
  int i;
  os << "ConminThetaVars" << endl; for(i=0;i<N1;i++) { os << conminThetaVars[i] << endl; }
  os << "ConminThetaLowerBnds" << endl; for(i=0;i<N1;i++) { os << conminThetaLowerBnds[i] << endl; }
  os << "ConminThetaUpperBnds" << endl; for(i=0;i<N1;i++) { os << conminThetaUpperBnds[i] << endl; }
  os << "ConstraintVector" << endl; for(i=0;i<conminSingleArray;i++) { os << constraintVector[i] << endl; }
  os << "SCAL" << endl; for(i=0;i<N1;i++) { os << SCAL[i] << endl; }
  os << "DF" << endl; for(i=0;i<N1;i++) { os << DF[i] << endl; }
  os << "A" << endl; for(i=0;i<N1*N3;i++) { os << A[i] << endl; }
  os << "S" << endl; for(i=0;i<N1;i++) { os << S[i] << endl; }
  os << "G1" << endl; for(i=0;i<N2;i++) { os << G1[i] << endl; }
  os << "G2" << endl; for(i=0;i<N2;i++) { os << G2[i] << endl; }
  os << "B" << endl; for(i=0;i<N3*N3;i++) { os << B[i] << endl; }
  os << "C" << endl; for(i=0;i<N4;i++) { os << C[i] << endl; }
  os << "ISC" << endl; for(i=0;i<N2;i++) { os << ISC[i] << endl; }
  os << "IC" << endl; for(i=0;i<N3;i++) { os << IC[i] << endl; }
  os << "MS1" << endl; for(i=0;i<N5;i++) { os << MS1[i] << endl; }
  os << "N1: " << N1 << endl;
  os << "N2: " << N2 << endl;
  os << "N3: " << N3 << endl;
  os << "N4: " << N4 << endl;
  os << "N5: " << N5 << endl;
  os << "N5: " << N5 << endl;
  os << "DELFUN: " << DELFUN << endl;
  os << "DABFUN: " << DABFUN << endl;
  os << "FDCH: " <<  FDCH<< endl;
  os << "FDCHM: " <<  FDCHM<< endl;
  os << "CT: " << CT << endl;
  os << "CTMIN: " << CTMIN << endl;
  os << "CTL: " << CTL << endl;
  os << "CTLMIN: " << CTLMIN << endl;
  os << "ALPHAX: " << ALPHAX << endl;
  os << "ABOBJ1: " << ABOBJ1 << endl;
  os << "THETA: " << THETA << endl;
  os << "maxLikelihoodEst: " << maxLikelihoodEst << endl;
  os << "dim: " << dim << endl;
  os << "numcon: " << numcon << endl;
  os << "NSIDE: " << NSIDE << endl;
  os << "IPRINT: " << IPRINT  << endl;
  os << "NFDG: " << NFDG  << endl;
  os << "NSCAL: " << NSCAL << endl;
  os << "LINOBJ: " << LINOBJ << endl;
  os << "ITMAX: " << ITMAX  << endl;
  os << "ITRM: " << ITRM << endl;
  os << "ICNDIR: " << ICNDIR  << endl;
  os << "IGOTO: " << IGOTO << endl;
  os << "NAC: " << NAC << endl;
  os << "conminInfo: " << conminInfo  << endl;
  os << "INFOG: " << INFOG  << endl;
  os << "ITER: " << ITER  << endl;
  os << "pts: " << pts  << endl;
  os << "numNewPts: " << numNewPts << endl;
  os << "xMatrix" << endl; for(i=0;i<pts;i++) { os << xMatrix[i] << endl; }
  os << "yValueVector" << endl; for(i=0;i<pts;i++) { os <<  yValueVector[i] << endl; }
  os << "xNewVector" << endl; for(i=0;i<pts;i++) { os <<  xNewVector[i]  << endl; }
  os << "yNewVector" << endl; for(i=0;i<pts;i++) { os <<  yNewVector[i] << endl; }
  os << "betaHat: " <<  betaHat << endl; 
  os << "rhsTermsVector" << endl; for(i=0;i<pts;i++) { os <<  rhsTermsVector[i] << endl; }
  os << "iPivotVector" << endl; for(i=0;i<pts;i++) { os <<  iPivotVector[i] << endl; }
  os << "correlationMatrix" << endl; for(i=0;i<pts*pts;i++) { os <<  correlationMatrix[i] << endl; }
  os << "invcorrelMatrix" << endl; for(i=0;i<pts*pts;i++) { os <<  invcorrelMatrix[i] << endl; }
  os << "fValueVector" << endl; for(i=0;i<pts;i++) { os <<  fValueVector[i] << endl; }
  os << "fRinvVector" << endl; for(i=0;i<pts;i++) { os <<  fRinvVector[i] << endl; }
  os << "yfbVector" << endl; for(i=0;i<pts;i++) { os <<  yfbVector[i] << endl; }
  os << "yfbRinvVector" << endl; for(i=0;i<pts;i++) { os <<  yfbRinvVector[i] << endl; }
  os << "rXhatVector" << endl; for(i=0;i<pts;i++) { os <<  rXhatVector[i] << endl; }
  os << "workVector" << endl; for(i=0;i<pts;i++) { os <<  workVector[i] << endl; }
  os << "workVectorQuad" << endl; for(i=0;i<4*pts;i++) { os <<  workVectorQuad[i] << endl; }
  os << "iworkVector" << endl; for(i=0;i<pts;i++) { os <<  iworkVector[i] << endl; }
  os << "numSampQuad: " << pts  << endl;
  os << "conminSingleArray: " << pts  << endl;
}

//_____________________________________________________________________________
// Helper methods 
//_____________________________________________________________________________

void KrigingSurface::buildModel()
{
  dim = static_cast<int>(dataItr->xSize());
  pts = static_cast<int>(dataItr->elementCount());
  unsigned i = 0;
  dataItr->toFront();
  while (!dataItr->isEnd()) {
    SurfPoint& current = dataItr->currentElement();
    for (int j = 0; j < dim; j++) {
       xMatrix[j*pts+i] = current.X()[j];
    }
    yValueVector[i] = current.F(dataItr->responseIndex());
    dataItr->nextElement();
    i++;
  } 
  if (runConminFlag) {
    // call to F77 driver code that runs CONMIN and the kriging software
    conminInfo = 0; // this flag must be zero before calling CONMIN
    IGOTO      = 0; // initialize CONMIN's internal flag
    iFlag      = 1; // flag for kriging code: value=1 used w/ optimization

    // Print out xMatrix
    //cout << "xMatrix before call to conmin" << endl;
    //for (int i = 0; i < dim*pts; i++) {
    //        cout << xMatrix[i] << endl;
    //}
    //cout << "yValueVector before call to conmin" << endl;
    //for (int j = 0; j < pts; j++) {
    //        cout << yValueVector[j] << endl;
    //}
#ifdef PRINT_DEBUG
    ofstream before("before.txt",ios::out);
    //ostream& os = cout;
    //os << "Before call to callconmin in buildModel" << endl;
    printConminVariables(before);
    before.close();
#endif
    callconmin(conminThetaVars[0], conminThetaLowerBnds[0],
	       conminThetaUpperBnds[0],
	       constraintVector[0],SCAL[0],DF[0],A[0],S[0],G1[0],G2[0],
	       B[0],C[0],ISC[0],IC[0],MS1[0],N1,N2,N3,N4,N5,
	       DELFUN,DABFUN,FDCH,FDCHM,CT,CTMIN,CTL,CTLMIN,ALPHAX,ABOBJ1,THETA,
	       maxLikelihoodEst,dim,numcon,NSIDE,IPRINT,NFDG,NSCAL,LINOBJ,
	       ITMAX,ITRM,ICNDIR,IGOTO,NAC,
	       conminInfo,INFOG,ITER,pts,numNewPts,iFlag, 
	       xMatrix[0],yValueVector[0],xNewVector[0],yNewVector[0],betaHat,
	       rhsTermsVector[0],iPivotVector[0],correlationMatrix[0],
	       invcorrelMatrix[0],fValueVector[0],fRinvVector[0], 
	       yfbVector[0],yfbRinvVector[0],rXhatVector[0],
	       workVector[0],workVectorQuad[0],iworkVector[0], 
	       numSampQuad,conminSingleArray);
#ifdef PRINT_DEBUG
    ofstream after("after.txt",ios::out);
    //os << "After call to callconmin in buildModel" << endl;
    printConminVariables(after);
    after.close();
#endif
    // copy CONMIN's theta variables to the thetaVector array for use in the
    // kriging model evaluation phase
    for (int i=0;i<dim;i++) {
      thetaVector[i] = conminThetaVars[i];
    }

    ////  run_conmin_flag = 0;
  }
  else {
    // Call to krigmodel() to bypass CONMIN. This approach uses the
    // user-supplied correlation values specified in the input file.
    iFlag = 1;
#ifdef PRINT_DEBUG
    ostream& os = cout;
    os << "Before call to krigmodel in buildModel" << endl;
    printKrigModelVariables(os);
#endif
    krigmodel(dim,pts,numNewPts,
	      iFlag,thetaVector[0],xMatrix[0],yValueVector[0],xNewVector[0],
	      yNewVector[0],betaHat,rhsTermsVector[0],maxLikelihoodEst,
	      iPivotVector[0],correlationMatrix[0],invcorrelMatrix[0],
	      fValueVector[0],fRinvVector[0],yfbVector[0],yfbRinvVector[0],
	      rXhatVector[0],workVector[0],workVectorQuad[0],iworkVector[0],
	      numSampQuad);
#ifdef PRINT_DEBUG
    os << "After call to krigmodel in buildModel" << endl;
    printKrigModelVariables(os);
#endif
  }
}

//_____________________________________________________________________________
// Testing 
//_____________________________________________________________________________

#ifdef __TESTING_MODE__
  int KrigingSurface::constructCount = 0;
  int KrigingSurface::copyCount = 0;
  int KrigingSurface::destructCount = 0;
#endif
