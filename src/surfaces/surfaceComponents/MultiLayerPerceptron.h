// These lines were automatically prepended
#ifndef SURFPACK_CONFIG_H
#define SURFPACK_CONFIG_H
#include "config.h"
#endif
// End of prepended lines

#ifndef MULTILAYERPERCEPTRON_H
#define MULTILAYERPERCEPTRON_H

#include "SurfPoint.h"
#include "SurfData.h"
#include "SurfDataIterator.h"
#include "SkipSurfDataIterator.h"
#include <vector>

class MultiLayerPerceptron 
{

public:
    MultiLayerPerceptron(SurfData* trainingData, int responseIndex, int numHiddenNodes, int maxTrainingIterations, double goalTrainingError);
    virtual ~MultiLayerPerceptron();
    virtual int getNumWeights() const;
    virtual void setWeights(const std::vector<double>& weightVector);
    virtual double onePattern(const std::vector<double>& x);
    virtual int getNumInputs() const;
    virtual int getNumHiddenNodes() const;
    virtual void load(std::string filename);
    virtual void save(std::string filename);
    virtual double oneEpoch(AbstractSurfDataIterator* iter = 0);
    virtual double trainNetwork() = 0;
    virtual void setTestData(SurfData* testData);
    
    virtual void detailedAnalysis();
    virtual double details(SurfData* data, std::ostream& os);
    
    static const int backprop;
    static const int particleSwarm;

protected:
//    int myrange;
//    int myprecision;
//    std::vector<double> eToX;

    SurfData* trainingData;
    int responseIndex;
    int numHiddenNodes;
    int maxTrainingIterations;
    double goalTrainingError;
    int numInputs;
    
    std::vector<double> inputs;
    double target;
    SurfData* testData;
    // data for network nodes by STL
    std::vector< std::vector<double> > weightsToHidden;
    std::vector< double > weightsToOutput;
    std::vector< double > hiddenActivations;
    double outputActivation;
    std::vector< std::vector<double> > deltaWeightsToHidden;
    std::vector< double > deltaWeightsToOutput;
    double deltaOutputNode;
    std::vector< double > deltaHiddenNodes;
    

    // helper functions for neural network

    void randomizeWeights();
    void initialize();
    double sigmoid(double net);
    double uniformRand();
    void fireNeurons();
    //void interpolate();
    void read(istream& is);
    void print(ostream& os);
    bool makeFrameThisIteration(int iteration);
    MultiLayerPerceptron() {}
};


#endif
