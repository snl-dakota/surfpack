// These lines were automatically prepended
#ifndef SURFPACK_CONFIG_H
#define SURFPACK_CONFIG_H
#include "config.h"
#endif
// End of prepended lines

#ifndef PARTICLESWARMMLP_H 
#define PARTICLESWARMMLP_H 

#include "MultiLayerPerceptron.h"
#include "OptimizationProblem.h"
#include "Swarm.h"
#include "Particle.h"

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
