/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#ifdef HAVE_CONFIG_H
#include "surfpack_config.h"
#endif

#include "BackpropMLP.h"
//#include "ANNWeights.h"
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

#define TRAIN_ON
//#define TEST_ON

BackpropMLP::BackpropMLP(SurfData* trainingData, int responseIndex, int numHiddenNodes, int maxTrainingIterations, double goalTrainingError, double learningRate)
    : MultiLayerPerceptron(trainingData, responseIndex, numHiddenNodes, maxTrainingIterations, goalTrainingError),
      learningRate(learningRate) 
{

}

BackpropMLP::~BackpropMLP()
{
}


double BackpropMLP::trainNetwork()
{

    int iterations = 0;
    double error = DBL_MAX;
    SurfDataIterator iter(trainingData);
    //while (iterations < maxTrainingIterations && error > goalTrainingError) {
    int maxIt = maxTrainingIterations;
    //double gError = goalTrainingError;
    while (iterations < maxIt ) { 
        //cout << "iterations: " << iterations << " maxIterations: " << maxTrainingIterations << endl;
	error = oneEpoch(&iter);
	int i = (iterations % 10 == 0) ? iterations / 10 : INT_MAX;
	if (makeFrameThisIteration(i)) {
		detailedAnalysis();
        }
	iterations++;
	if (iterations % 1000 == 0) cout << setw(8) << iterations << setw(15) << error << endl;
    }
    return error;
}

double BackpropMLP::oneEpoch(AbstractSurfDataIterator* iter)
{
	bool ownIter = false;
	if (!iter) {
	    iter = new SurfDataIterator(trainingData);
	    ownIter = true;
	}
	double sse = 0.0;  // sum squared error
	iter->shuffle();
	SurfPoint* current;
	while(!iter->isEnd()) {
		current = iter->currentElement();
		inputs = current->getX();
		target = current->getF(responseIndex);
		fireNeurons();
		backpropagateErrors();
		sse += (target - outputActivation) * (target - outputActivation);
		iter->nextElement();
	}
	if (ownIter) {
	    delete iter;
	}
	return 0.5 * sse; // the .5 makes the bp-rule math come out nice
}



void BackpropMLP::backpropagateErrors()
{
    int hiddenNode, inputNode; // counters;
    deltaOutputNode =  (target - outputActivation);
    for (hiddenNode = 0; hiddenNode < numHiddenNodes; hiddenNode++) {
        deltaHiddenNodes[hiddenNode] = hiddenActivations[hiddenNode] * (1.0 - hiddenActivations[hiddenNode]) * 
	    (weightsToOutput[hiddenNode] * deltaOutputNode);
        weightsToOutput[hiddenNode] += learningRate * deltaOutputNode * hiddenActivations[hiddenNode];	
    }
    weightsToOutput[numHiddenNodes] += learningRate * deltaOutputNode;
    
    // update weights from input to hidden layer
    for (hiddenNode = 0; hiddenNode < numHiddenNodes; hiddenNode++) {
	for (inputNode = 0; inputNode < numInputs; inputNode++) {
	    weightsToHidden[hiddenNode][inputNode] += learningRate * deltaHiddenNodes[hiddenNode] * inputs[inputNode];
	}
        weightsToHidden[hiddenNode][numInputs] += learningRate * deltaHiddenNodes[hiddenNode];
    }
}    


