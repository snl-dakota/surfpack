// ----------------------------------------------------------
// Project: SURFPACK++
//
// File:        KrigingSurface.h
// Author:      Mark Richards	 
// Created:     June 3, 2003
// Modified:    
//
// Description: 
// + The KrigingSurface class does Kriging interpolation 
// (similar to splines) to fit on the points
// ----------------------------------------------------------

#ifndef __KRIGING_SURFACE_H__
#define __KRIGING_SURFACE_H__

#include <iostream>
#include <vector>
#include <string>

#include "SurfException.h"
#include "SurfPoint.h"
#include "Surface.h"
#include "SkipSurfDataIterator.h"
#include "SurfDataIterator.h"


class KrigingSurface : public Surface
{

public:

// constructor/destructor
////////////////////////////////

   //KrigingSurface(const std::vector<double> &);
   KrigingSurface(SurfData* surfData, int responseIndex = 0);
   KrigingSurface(std::istream& is);
   KrigingSurface(std::string filename);
//   KrigingSurface(SurfData *,const std::vector<double> &);
//   KrigingSurface(const KrigingSurface &);
   
   // a do nothing destructor
   //
   virtual ~KrigingSurface(); 

// member functions
////////////////////////////////


   virtual int getMinPointCount(int dim) const;
   virtual int getDimension() const;
   virtual std::ostream& write(std::ostream & os = std::cout); 
   virtual std::string getType() const;
   virtual void save(std::string filename);

   /// Save the surface to a file in binary format
   virtual void saveBinary(std::string filename);

   /// Load the surface from a file
   virtual void loadBinary(std::string filename);

private: 
   KrigingSurface();
   std::ostream& printKrigEvalVars(ostream& os);
   std::ostream& printKrigModelVariables(std::ostream& os);
   std::ostream& printConminVariables(std::ostream& os);

protected:

   bool needsCleanup;
   int numdv;
   int numsamp;
   int dim;
   int pts;
   virtual double calculate(const std::vector<double> & x) throw(SurfException);
   virtual void calculateInternals(AbstractSurfDataIterator* iter);
   void initialize();
   void cleanup();
   void usePreComputedCorrelationVector(std::vector<double> vals) throw(SurfException);
   virtual void build(AbstractSurfDataIterator* iter) throw(SurfException);
   void buildModel();
   bool runConminFlag;



    //
    //- Heading: Data
    //

  /*******************CONMIN DATA**************************/

  /// Size variable for CONMIN arrays. See CONMIN manual.
  /** N1 = number of variables + 2 */
  int N1;
  /// Size variable for CONMIN arrays. See CONMIN manual.
  /** N2 = number of constraints + 2*(number of variables) */
  int N2;
  /// Size variable for CONMIN arrays. See CONMIN manual.
  /** N3 = Maximum possible number of active constraints.*/
  int N3;
  /// Size variable for CONMIN arrays. See CONMIN manual.
  /** N4 = Maximum(N3,number of variables) */
  int N4;
  /// Size variable for CONMIN arrays. See CONMIN manual.
  /** N5 = 2*(N4) */
  int N5;
  /// Array size parameter needed in interface to CONMIN.
  int conminSingleArray;
  /// CONMIN variable: Number of constraints.
  int numcon;
  /// CONMIN variable: Finite difference flag.
  int NFDG;
  /// CONMIN variable: Flag to control amount of output data.
  int IPRINT;  
  /// CONMIN variable: Flag to specify the maximum number of iterations.
  int ITMAX;  
  /// CONMIN variable: Relative finite difference step size.
  double FDCH; 
  /// CONMIN variable: Absolute finite difference step size.
  double FDCHM;
  /// CONMIN variable: Constraint thickness parameter 
  /** The value of CT decreases in magnitude during optimization.*/
  double CT;
  /// CONMIN variable: Minimum absolute value of CT used during optimization.
  double CTMIN;
  /// CONMIN variable: Constraint thickness parameter for linear and
  /// side constraints.
  double CTL;
  /// CONMIN variable: Minimum value of CTL used during optimization.
  double CTLMIN;
  /// CONMIN variable: Relative convergence criterion threshold.
  /*** Threshold for the minimum relative change in the objective function. */
  double DELFUN; 
  /// CONMIN variable: Absolute convergence criterion threshold.
  /*** Threshold for the minimum relative change in the objective function. */
  double DABFUN;
  /// CONMIN variable: status flag for optimization
  int conminInfo;
  /// Internal CONMIN array.
  /** Move direction in N-dimensional space.*/
  double *S;
  /// Internal CONMIN array.
  /** Temporary storage of constraint values.*/
  double *G1;
  /// Internal CONMIN array.
  /** Temporary storage of constraint values.*/
  double *G2;
  /// Internal CONMIN array.
  /** Temporary storage for computations involving array S.*/
  double *B;
  /// Internal CONMIN array.
  /** Temporary storage for use with arrays B and S.*/
  double *C;
  /// Internal CONMIN array.
  /** Temporary storage for use with arrays B and S.*/
  int    *MS1;
  /// Internal CONMIN array.
  /** Vector of scaling parameters for design parameter values.*/
  double *SCAL;
  /// Internal CONMIN array.
  /** Temporary storage for analytic gradient data.*/
  double *DF;
  /// Internal CONMIN array.
  /** Temporary 2-D array for storage of constraint gradients.*/
  double *A;
  /// Internal CONMIN array.
  /** Array of flags to identify linear constraints. (not used in this
      implementation of CONMIN) */
  int    *ISC; 
  /// Internal CONMIN array.
  /** Array of flags to identify active and violated constraints */
  int    *IC;

  /// Temporary array of design variables used by CONMIN (length N1 = numdv+2)
  double *conminThetaVars;
  /// Temporary array of lower bounds used by CONMIN (length N1 = numdv+2)
  double *conminThetaLowerBnds;
  /// Temporary array of upper bounds used by CONMIN (length N1 = numdv+2)
  double *conminThetaUpperBnds;

  /// Internal CONMIN variable: 1-D search parameter.
  double ALPHAX;
  /// Internal CONMIN variable: 1-D search parameter.
  double ABOBJ1;
  /// Internal CONMIN variable: mean value of push-off factor.
  double THETA;
  /// Internal CONMIN variable: "participation coefficient".
  double PHI;
  /// Internal CONMIN variable: side constraints parameter.
  int  NSIDE;
  /// Internal CONMIN variable: scaling control parameter.
  int  NSCAL;
  /// Internal CONMIN variable: estimate of 1+(max # of active constraints).
  int  NACMX1;
  /// Internal CONMIN variable: linear objective function identifier (unused).
  int  LINOBJ;
  /// Internal CONMIN variable: diminishing return criterion iteration number.
  int  ITRM;
  /// Internal CONMIN variable: conjugate direction restart parameter.
  int  ICNDIR;
  /// Internal CONMIN variable: internal optimization termination flag.
  int  IGOTO;
  /// Internal CONMIN variable: number of active and violated constraints.
  int  NAC;
  /// Internal CONMIN variable: gradient information flag.
  int  INFOG;
  /// Internal CONMIN variable: iteration count.
  int  ITER;


  /********************KRIGING CODE DATA****************************/
  
  /// Fortran77 flag for kriging computations.
  /** iFlag=1 computes vector and matrix terms for the kriging surface, 
      iFlag=2 computes the response value (using kriging) at the 
      user-supplied design point.*/
  int  iFlag;
  /// Estimate of the beta term in the kriging model..
  double betaHat;
  /// Error term computed via Maximum Likelihood Estimation.
  double maxLikelihoodEst;
  /// Size variable for the arrays used in kriging computations.
  int numNewPts;
  /// Size variable for the arrays used in kriging computations.
  int numSampQuad;
  /// Array of correlation parameters for the kriging model.
  double *thetaVector;
  /// A 2-D array of design points used to build the kriging model.
  double *xMatrix;
  /// A vector representing a data point.
  double *xArray;
  /// Array of response values corresponding to the array of design points.
  double *yValueVector;
  /// A 2-D array of design points where the kriging model will be evaluated.
  double *xNewVector;
  /// Array of response values corresponding to the design points
  /// specified in xNewVector.
  double *yNewVector;
  /// Array of lower bounds in optimizer-to-kriging interface.
  double *thetaLoBndVector;
  /// Array of upper bounds in optimizer-to-kriging interface.
  double *thetaUpBndVector;
  /// Array of constraint values (used with optimizer)
  double *constraintVector;
  /// Internal array for kriging Fortran77 code: matrix algebra result.
  double *rhsTermsVector;
  /// Internal array for kriging Fortran77 code: pivot vector for
  /// linear algebra.
  int  *iPivotVector;
  /// Internal array for kriging Fortran77 code: correlation matrix.
  double *correlationMatrix;
  /// Internal array for kriging Fortran77 code: inverse correlation matrix.
  double *invcorrelMatrix;
  /// Internal array for kriging Fortran77 code: response value vector.
  double *fValueVector;
  /// Internal array for kriging Fortran77 code: vector*matrix result.
  double *fRinvVector;
  /// Internal array for kriging Fortran77 code: vector arithmetic result.
  double *yfbVector;
  /// Internal array for kriging Fortran77 code: vector*matrix result.
  double *yfbRinvVector;
  /// Internal array for kriging Fortran77 code: local correlation vector.
  double *rXhatVector;
  /// Internal array for kriging Fortran77 code: temporary storage.
  double *workVector;
  /// Internal array for kriging Fortran77 code: temporary storage.
  double *workVectorQuad;
  /// Internal array for kriging Fortran77 code: temporary storage.
  int  *iworkVector;
};

#endif
