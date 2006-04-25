/*  _______________________________________________________________________

    DAKOTA: Design Analysis Kit for Optimization and Terascale Applications
    Copyright (c) 2001, Sandia National Laboratories.

    Surfpack: A Software Library of Multidimensional Surface Fitting Methods

    Surfpack is distributed under the DAKOTA GNU General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */
#include "surfpack_config.h"

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
