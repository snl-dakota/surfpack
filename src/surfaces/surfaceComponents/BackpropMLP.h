// These lines were automatically prepended
#ifndef SURFPACK_CONFIG_H
#define SURFPACK_CONFIG_H
#include "config.h"
#endif
// End of prepended lines

#ifndef BACKPROPMLP_H 
#define BACKPROPMLP_H 

#include "MultiLayerPerceptron.h"

#include "SurfPoint.h"
#include "SurfData.h"
#include <vector>

class BackpropMLP : public MultiLayerPerceptron
{

public:
    BackpropMLP(SurfData* trainingData, int responseIndex, int numHiddenNodes, int maxTrainingIterations, double goalTrainingError, double learningRate = .01);
    virtual ~BackpropMLP();
    virtual double trainNetwork();
    virtual double oneEpoch(AbstractSurfDataIterator* iter);
    virtual void backpropagateErrors();

protected:
    double learningRate;
    
private:
    BackpropMLP() {}
};


#endif
