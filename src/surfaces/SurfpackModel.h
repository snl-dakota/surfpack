/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#ifndef __SURFPACK_MODEL_H__
#define __SURFPACK_MODEL_H__
#ifdef HAVE_CONFIG_H
#include "surfpack_config.h"
#endif

#include "surfpack_system_headers.h"
#include "ModelScaler.h"

class SurfData;

///////////////////////////////////////////////////////////
///	Surfpack Model Parameters 
///////////////////////////////////////////////////////////

typedef std::pair< std::string, std::string > ModelParam;
typedef std::map< std::string, std::string> ParamMap;

///////////////////////////////////////////////////////////
///	Surfpack Model 
///////////////////////////////////////////////////////////

class SurfpackModel
{
public:
  SurfpackModel(unsigned ndims_in);
  SurfpackModel(const SurfpackModel& other);
  virtual VecDbl operator()(const SurfData& data) const;
  double operator()(const VecDbl& x) const;
  virtual VecDbl gradient(const VecDbl& x) const = 0;
  virtual std::string asString() const = 0;
  virtual ~SurfpackModel();
  virtual void scaler(ModelScaler* ms);
  ModelScaler* scaler() const;
  unsigned size() const { return ndims;}
protected:
  virtual double evaluate(const VecDbl& x) const = 0;
  unsigned ndims;
  ModelScaler*  mScaler;
};

///////////////////////////////////////////////////////////
///	Surfpack Model Factory
///////////////////////////////////////////////////////////

class SurfpackModelFactory
{

public:
  SurfpackModelFactory();
  SurfpackModelFactory(const ParamMap& args);
  virtual SurfpackModel* Create(const SurfData& sd) = 0;
  virtual SurfpackModel* Create(const std::string& model_string) = 0;
  virtual unsigned minPointsRequired();
  virtual void config();
  const ParamMap& parameters() const;
  void add(std::string name, std::string value);
protected:
  ParamMap params;
  unsigned ndims;
};

#endif
