#include "ParticleSwarmMLP.h"
#include "ANNWeights.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <cfloat>
#include <ctime>
#include <algorithm>
using namespace std;

ParticleSwarmMLP::ParticleSwarmMLP(SurfData* trainingData, int responseIndex, int numHiddenNodes, int maxTrainingIterations, double goalTrainingError, int swarmSize, int sociometry)
    : MultiLayerPerceptron(trainingData, responseIndex, numHiddenNodes, maxTrainingIterations, goalTrainingError),
      swarmSize(swarmSize), sociometry(sociometry)

{
    weightsOpt = new ANNWeights(this);
}

ParticleSwarmMLP::~ParticleSwarmMLP()
{
    delete weightsOpt;
}


double ParticleSwarmMLP::trainNetwork()
{
    int trainingIntervals = maxTrainingIterations / swarmSize;
    if (maxTrainingIterations % swarmSize != 0) {
	    cerr << "Wasted function evals" << endl;
    }
    Swarm* swarm = new Swarm(weightsOpt, swarmSize, sociometry);
    vector<double> bestWeights(getNumWeights());
    vector< vector<double> > dummy;
    double error = 0.0;
    for (int i = 0; i < trainingIntervals; i++) {
        error = swarm->optimize(bestWeights, swarmSize, goalTrainingError, dummy);
        setWeights(bestWeights);
        if (makeFrameThisIteration(i)) {
                cout << i*swarmSize << "   ";
                detailedAnalysis();
        }
    if (i*swarmSize % 1000 == 0) cout << setw(10) << i*swarmSize << setw(15) << error << endl;
    }
    //cout << setw(10) << maxTrainingIterations << setw(15) << error << endl;
    delete swarm;
    return error; 
}

