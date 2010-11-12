#ifndef _SURFPACKMODEL_H_
#define _SURFPACKMODEL_H_

#include "NKM_SurfMat.hpp"
#include "NKM_Optimize.hpp"
#include <iostream>
#include <exception>

namespace nkm {

using std::cerr;
using std::endl;

class SurfPackModel 
{
protected:
  SurfData sdBuild;  
  SurfDataScaler scaler;

public:
  
  SurfPackModel() : sdBuild(), scaler(sdBuild) {};

  SurfPackModel(const SurfData& sd,int jout_keep) : sdBuild(sd,jout_keep), scaler(sdBuild) {};

  /// adjust correlations to be feasible with respect to condition
  /// number constraints
  virtual MtxDbl& makeGuessFeasible(MtxDbl& correlations, 
				    OptimizationProblem *opt)
  {
    // default at base class is no-op
  }

  virtual void getRandGuess(MtxDbl& guess) const{


  }

  virtual void set_conmin_parameters(OptimizationProblem& opt) const{
  };

  virtual void set_direct_parameters(OptimizationProblem& opt) const{
  };
  

  /// the objective function, i.e. the negative log(likelihood);
  /// minimizing this produces a "good" KrigingModel)
  virtual double objective(const MtxDbl& correlations)
  {
    cerr << "Derived class does not implement objective" << endl;
    throw(std::string("Derived does not implement"));
  }

  /// the objective function, i.e. the negative log(likelihood), and
  /// its gradient; minimizing the objective function produces a good
  /// KrigingModel
  virtual void objectiveAndGradient(double& Obj, MtxDbl& GradObj,
				    const MtxDbl& correlations)
  {
    cerr << "Derived class does not implement objectiveAndGradient" << endl;
    throw(std::string("Derived does not implement"));
  }  


  /// objective plus condition number constraints
  virtual void objectiveAndConstraints(double& Obj, MtxDbl& Con, 
				       const MtxDbl& correlations)
  {
    cerr << "Derived class does not implement objectiveAndConstraints" << endl;
    throw(std::string("Derived does not implement"));
  }


  /// objective plus condition number constraints with gradients
  virtual void objectiveAndConstraintsAndGradients(double& Obj, MtxDbl& Con, 
						   MtxDbl& GradObj, 
						   MtxDbl& GradCon, 
						   const MtxDbl& correlations)
  {
    cerr << "Derived class does not implement objectiveAndConstraintsAndGradients" << endl;
    throw(std::string("Derived does not implement"));
  }



};

} // end namespace nkm

#endif
