/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
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
