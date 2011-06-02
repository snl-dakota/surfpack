#ifndef _SURFPACKMODEL_H_
#define _SURFPACKMODEL_H_

#include "NKM_SurfMat.hpp"
#include "NKM_Optimize.hpp"
#include <iostream>
#include <exception>

namespace nkm {

class SurfPackModel 
{
protected:
  SurfData sdBuild;  
  SurfDataScaler scaler;

public:
  
  SurfPackModel() : sdBuild(), scaler(sdBuild) {};

  SurfPackModel(const SurfData& sd,int jout_keep) : sdBuild(sd,jout_keep), scaler(sdBuild) {};

  virtual void create() {
    std::cerr << "the create() function has not been implemented for this model type" << std::endl;
    return;
  };

  virtual std::string model_summary_string() const {
    std::string mod_sum_str="the model_summary_string() function has not been implemented for this model\n";
    return mod_sum_str;
  };


  virtual double evaluate(const  MtxDbl& xr) const =0;

  virtual MtxDbl& evaluate(MtxDbl& y, const MtxDbl& xr) const
  {
    int nrowsxr=xr.getNRows();
    int ncolsxr=xr.getNCols();
    assert((ncolsxr==sdBuild.getNVarsr())&&(nrowsxr>0));
    y.newSize(nrowsxr,1);

    if(nrowsxr==1) {
      y(0,0)=evaluate(xr);
      return y;
    }
      
    MtxDbl xr_temp(1,ncolsxr);
    for(int ipt=0; ipt<nrowsxr; ++ipt) {
      xr.getRows(xr_temp,ipt);
      y(ipt,0)=evaluate(xr_temp);
    }
    return y;
  };

  virtual double eval_variance(const MtxDbl& xr) const {
    std::cerr << "This model doesn't have an implemented function to return a variance" << std::endl;
    assert(false);
    return (0.0/0.0);
  };

  virtual MtxDbl& eval_variance(MtxDbl& var, const MtxDbl& xr) const
  {
    int nrowsxr=xr.getNRows();
    int ncolsxr=xr.getNCols();
    assert((ncolsxr==sdBuild.getNVarsr())&&(nrowsxr>0));
    var.newSize(nrowsxr,1);

    if(nrowsxr==1) {
      var(0,0)=eval_variance(xr);
      return var;
    }
      
    MtxDbl xr_temp(1,ncolsxr);
    for(int ipt=0; ipt<nrowsxr; ++ipt) {
      xr.getRows(xr_temp,ipt);
      var(ipt,0)=eval_variance(xr_temp);
    }
    return var;
  };

  virtual MtxDbl& evaluate_d1y(MtxDbl& d1y, const MtxDbl& xr) const =0;

  virtual MtxDbl& evaluate_d2y(MtxDbl& d2y, const MtxDbl& xr) const =0;


  /// adjust correlations to be feasible with respect to condition
  /// number constraints
  virtual MtxDbl& makeGuessFeasible(MtxDbl& correlations, 
				    OptimizationProblem *opt)
  {
    // default at base class is no-op
    return correlations;
  };

  virtual void getRandGuess(MtxDbl& guess) const{


  };

  virtual void set_conmin_parameters(OptimizationProblem& opt) const{
  };

  virtual void set_direct_parameters(OptimizationProblem& opt) const{
  };
  

  /// the objective function, i.e. the negative log(likelihood);
  /// minimizing this produces a "good" KrigingModel)
  virtual double objective(const MtxDbl& correlations)
  {
    std::cerr << "Derived class does not implement objective" << std::endl;
    throw(std::string("Derived does not implement"));
  }

  /// the objective function, i.e. the negative log(likelihood), and
  /// its gradient; minimizing the objective function produces a good
  /// KrigingModel
  virtual void objectiveAndGradient(double& Obj, MtxDbl& GradObj,
				    const MtxDbl& correlations)
  {
    std::cerr << "Derived class does not implement objectiveAndGradient" << std::endl;
    throw(std::string("Derived does not implement"));
  }  


  /// objective plus condition number constraints
  virtual void objectiveAndConstraints(double& Obj, MtxDbl& Con, 
				       const MtxDbl& correlations)
  {
    std::cerr << "Derived class does not implement objectiveAndConstraints" << std::endl;
    throw(std::string("Derived does not implement"));
  }


  /// objective plus condition number constraints with gradients
  virtual void objectiveAndConstraintsAndGradients(double& Obj, MtxDbl& Con, 
						   MtxDbl& GradObj, 
						   MtxDbl& GradCon, 
						   const MtxDbl& correlations)
  {
    std::cerr << "Derived class does not implement objectiveAndConstraintsAndGradients" << std::endl;
    throw(std::string("Derived does not implement"));
  }



};

} // end namespace nkm

#endif
