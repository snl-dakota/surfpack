/*  _______________________________________________________________________

    DAKOTA: Design Analysis Kit for Optimization and Terascale Applications
    Copyright (c) 2001, Sandia National Laboratories.

    Surfpack: A Software Library of Multidimensional Surface Fitting Methods

    Surfpack is distributed under the DAKOTA GNU General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */
#include "surfpack_config.h"

#ifndef PARTICLESWARMMLP_H 
#define PARTICLESWARMMLP_H 

#include "MultiLayerPerceptron.h"
//#include "OptimizationProblem.h"
//#include "Swarm.h"
//#include "Particle.h"

#include "SurfPoint.h"
#include "SurfData.h"
#include <vector>

class ANNWeights;

class ParticleSwarmMLP : public MultiLayerPerceptron
{

public:
    ParticleSwarmMLP(SurfData* trainingData, int responseIndex, int numHiddenNodes, int maxTrainingIterations, double goalTrainingError, int swarmSize = 10, int sociometry = 2);
    virtual ~ParticleSwarmMLP();
    virtual double trainNetwork();

protected:
    int swarmSize;
    int sociometry;
    ANNWeights* weightsOpt;
private:
    ParticleSwarmMLP() {}
};


#endif
