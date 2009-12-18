/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#ifndef __CONMIN_H__ 
#define __CONMIN_H__ 

#ifdef HAVE_CONFIG_H
#include "surfpack_config.h"
#endif
#include "surfpack_system_headers.h"

#include "SurfpackModel.h"

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
             int& numdv, int& ncon, int& nside, int& iprint, int& nfdg,
             int& nscal, int& linobj, int& itmax, int& itrm, int& incdir,
             int& igoto, int& nac, int& info, int& infog, int& iter);




class Conmin {
public:
  Conmin(unsigned ndv_in);
  void bounds(const VecDbl& lower_bounds, const VecDbl& upper_bounds);
  virtual void optimize(VecDbl& x, double& final_val, unsigned max_iter) = 0;
  virtual double objective(const VecDbl& x) = 0;
  virtual VecDbl gradient(const VecDbl& x) = 0;
  virtual ~Conmin();
protected:
  VecDbl upperBounds;
  VecDbl lowerBounds;
  int NSIDE;
  unsigned ndv;
  
};

#endif
