/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
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
