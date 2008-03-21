/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#ifndef __MARS_MODEL_H__
#define __MARS_MODEL_H__
#ifdef HAVE_CONFIG_H
#include "surfpack_config.h"
#endif

#include "surfpack_system_headers.h"
#include "SurfpackModel.h"
#include "SurfpackMatrix.h"
class SurfPoint;

typedef float real;

#define MARS_F77 F77_FUNC(mars,MARS)
#ifdef __cplusplus
extern "C" /* prevent C++ name mangling */
#endif
void MARS_F77(int&, int&, real*, real*, real*, int&, int&, int*,
	      real*, int*, real*, double*, int*);

#define FMODM_F77 F77_FUNC(fmodm,FMODM)
#ifdef __cplusplus
extern "C" /* prevent C++ name mangling */
#endif
void FMODM_F77(int&, int&, real*, real*, int*, real*, real*);

class MarsModel : public SurfpackModel
{
public:
  MarsModel(const unsigned dims, real* fm_in, int fmsize, int* im_in, 
    int imsize, int interp);
  virtual VecDbl gradient(const VecDbl& x) const;
  virtual std::string asString() const;
protected:
  virtual double evaluate(const VecDbl& x) const;
  VecDbl coeffs;

  std::vector<real> fm;
  std::vector<int> im;
  int interpolation;
friend class MarsModelTest;
};

///////////////////////////////////////////////////////////
///	Linear Regression Model Factory	
///////////////////////////////////////////////////////////

class MarsModelFactory : public SurfpackModelFactory 
{

public:
  MarsModelFactory();
  MarsModelFactory(const ParamMap& args);
  virtual SurfpackModel* Create(const SurfData& sd);
  virtual SurfpackModel* Create(const std::string& model_string);
  virtual void config();
protected:
  real* xMatrix;
  real* fm;
  int* im;
  int n;
  int np;
  int max_bases; 
  int max_interactions;
  int interpolation;
};
#endif
