// These lines were automatically prepended
#ifndef SURFPACK_CONFIG_H
#define SURFPACK_CONFIG_H
#include "config.h"
#endif
// End of prepended lines

#include "MultiLayerPerceptron.h"
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


const double PI = 2*acos(0);

const int MultiLayerPerceptron::backprop = 1;
const int MultiLayerPerceptron::particleSwarm = 2;

MultiLayerPerceptron::MultiLayerPerceptron(SurfData* trainingData, int responseIndex, int numHiddenNodes, int maxTrainingIterations, double goalTrainingError)
    : trainingData(trainingData), responseIndex(responseIndex), numHiddenNodes(numHiddenNodes), maxTrainingIterations(maxTrainingIterations), goalTrainingError(goalTrainingError)
{
    if (!trainingData) {
	    cerr << "No training patterns" << endl;
    } else if (trainingData->getDimension() <= 0) {
	   cerr << "Bad dimension of dataset in MLP::MLP" << endl;
    } else {
	    numInputs = trainingData->getDimension();
    }
	   
    if (numHiddenNodes < 1) {
	    this->numHiddenNodes = 10;
    } 
    
    initialize();
    randomizeWeights();
}

MultiLayerPerceptron::~MultiLayerPerceptron()
{
}

int MultiLayerPerceptron::getNumInputs() const
{
	return numInputs;
}

int MultiLayerPerceptron::getNumHiddenNodes() const
{
	return numHiddenNodes;
}


int MultiLayerPerceptron::getNumWeights() const
{
    return numHiddenNodes * (numInputs + 1) + numHiddenNodes + 1;
}

double MultiLayerPerceptron::oneEpoch(AbstractSurfDataIterator* iter)
{
	bool ownIter = false;
	if (!iter) {
	    iter = new SurfDataIterator(trainingData);
	    ownIter = true;
	}
	double sse = 0.0;  // sum squared error
	iter->toFront();
	SurfPoint* current;
	while(!iter->isEnd()) {
		current = iter->currentElement();
		inputs = current->getX();
		target = current->getF(responseIndex);
		fireNeurons();
		sse += (target - outputActivation) * (target - outputActivation);
		iter->nextElement();
	}
	if (ownIter) {
	    delete iter;
	}
	return 0.5 * sse; // the .5 makes the bp-rule math come out nice
}

double MultiLayerPerceptron::onePattern(const vector<double>& x)
{
    inputs = x;
    fireNeurons();
    return outputActivation;
}

void MultiLayerPerceptron::setWeights(const vector<double>& weightVector) 
{
    if (static_cast<unsigned>(getNumWeights()) != weightVector.size()) {
	    cerr << "Mismatch in MultiLayerPerceptron::setWeights" << endl;
    }
    else
    {
        int index = 0;
	int hiddenNode, inputNode; // counters
	for (hiddenNode = 0; hiddenNode < numHiddenNodes; hiddenNode++) {
		for (inputNode = 0; inputNode < numInputs; inputNode++) {
			weightsToHidden[hiddenNode][inputNode] = weightVector[index++];
		}
		weightsToHidden[hiddenNode][numInputs] = weightVector[index++]; // bias term
		weightsToOutput[hiddenNode] = weightVector[index++];
	}
	weightsToOutput[numHiddenNodes] = weightVector[index++]; // bias term
	if (static_cast<unsigned>(index) != weightVector.size()) {
		cerr << "Didn't count weights correctly in setWeights" << endl;
	}
    }
}

double MultiLayerPerceptron::sigmoid(double net)
{
    //return tanh(net);
	//if (net < -(myrange *.5)) {
	//	return -1.0;
	//} else if (net > (myrange*.5)) {
	//	return 1.0;
	//} else {
	//	int value = static_cast<int>(net*myprecision) + myrange*myprecision;
	//	return eToX[value];		
	//}
	//static const double c = 1.0 / (0.4 * sqrt(2.0*PI));
	//return c * exp(-.5*net*net);	
	return 1.0 / (1.0 + exp(-net));

	
}

double MultiLayerPerceptron::uniformRand()
{
    return static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
}

void MultiLayerPerceptron::initialize()
{
    weightsToHidden.resize(numHiddenNodes);
    for(vector< vector<double> >::iterator it = weightsToHidden.begin();
	it != weightsToHidden.end();
	++it ) {
	    it->resize(numInputs + 1);  // include input from bias term
    }
    weightsToOutput.resize(numHiddenNodes + 1); // include a weight from the hidden layer bias term
    hiddenActivations.resize(numHiddenNodes);
    deltaWeightsToHidden.resize(numHiddenNodes);
    for(vector< vector<double> >::iterator it2 = deltaWeightsToHidden.begin();
	it2 != deltaWeightsToHidden.end();
	++it2 ) {
	    it2->resize(numInputs + 1);  
    }
    deltaWeightsToOutput.resize(numHiddenNodes+1);
    deltaHiddenNodes.resize(numHiddenNodes);
}

void MultiLayerPerceptron::fireNeurons()
{
    int inputNode, hiddenNode; // counters
    //const vector<double>& inputs = currentPattern->getX();
    for(hiddenNode = 0; hiddenNode < numHiddenNodes; hiddenNode++) {
	hiddenActivations[hiddenNode] = 0.0;
	for(inputNode = 0; inputNode < numInputs; inputNode++) {
	    hiddenActivations[hiddenNode] += inputs[inputNode]*weightsToHidden[hiddenNode][inputNode];
	}
	hiddenActivations[hiddenNode] += weightsToHidden[hiddenNode][numInputs]; // add in the bias term
	hiddenActivations[hiddenNode] = sigmoid(hiddenActivations[hiddenNode]);
    }
    outputActivation = 0;
    for(hiddenNode = 0; hiddenNode < numHiddenNodes; hiddenNode++) {
	outputActivation += hiddenActivations[hiddenNode]*weightsToOutput[hiddenNode];
    }
    outputActivation += weightsToOutput[numHiddenNodes]; // add in the bias output
}



void MultiLayerPerceptron::randomizeWeights()
{
    int inputNode, hiddenNode; // counters
    
    //myrange = 40;
    //myprecision = 10000;
    //eToX.resize(myrange*myprecision);
    //double value;
    //for (unsigned i = 0; i < eToX.size(); i++) { 
    //    value = (static_cast<double>(i) / static_cast<double>(myprecision)) - .5 * static_cast<double>(myrange) ;
    //    eToX[i] = 1.0 / (1.0 + exp(-value));
    //    //cout << setw(6) << i << setw(10) << value << setw(15) << eToX[i] << endl;
    //}
   
    for(hiddenNode = 0; hiddenNode < numHiddenNodes; hiddenNode++) {
	for(inputNode = 0; inputNode < numInputs; inputNode++) {
		weightsToHidden[hiddenNode][inputNode] =0.1 * (rand() / static_cast<double>(RAND_MAX)) - .05;
	}
	weightsToHidden[hiddenNode][numInputs] = 0.1 * (rand() / static_cast<double>(RAND_MAX)) - .05; 
	weightsToOutput[hiddenNode] = 0.1 * (rand() / static_cast<double>(RAND_MAX)) - .05;
    }
    weightsToOutput[numHiddenNodes] =0.1 * (rand() / static_cast<double>(RAND_MAX)) - .05; 
}

void MultiLayerPerceptron::read(istream& is)
{
	int inputNode, hiddenNode; // counters
	for (hiddenNode = 0; hiddenNode < numHiddenNodes; hiddenNode++) {
		for (inputNode = 0; inputNode <= numInputs; inputNode++) {
			is >> weightsToHidden[hiddenNode][inputNode] ; 
		}
	}
        for (hiddenNode = 0; hiddenNode <= numHiddenNodes; hiddenNode++) {
	    is >> weightsToOutput[hiddenNode] ; 
	}
}

void MultiLayerPerceptron::print(ostream& os)
{
	int inputNode, hiddenNode; // counters
    
	for (hiddenNode = 0; hiddenNode < numHiddenNodes; hiddenNode++) {
		for (inputNode = 0; inputNode <= numInputs; inputNode++) {
			os << setprecision(20) << weightsToHidden[hiddenNode][inputNode] << endl; 
		}
	}
        for (hiddenNode = 0; hiddenNode <= numHiddenNodes; hiddenNode++) {
	    os << setprecision(20) << weightsToOutput[hiddenNode] << endl; 
	}
	    	
}


void MultiLayerPerceptron::load(string filename)
{
   ifstream infile("Weights/weights.txt",ios::in);
   read(infile);
   infile.close();
}

void MultiLayerPerceptron::save(string filename)
{
   ofstream outfile("Weights/weights.txt",ios::out);
   print(outfile);
   outfile.close();
}
   
void MultiLayerPerceptron::detailedAnalysis()
{
    static int called = 0;
    //return; // to eliminate the verbose output
    
    if (!testData) {
	    cerr << "No test Data in detailedAnalysis" << endl;
    }
    ostringstream fnstream ;
    fnstream << "analysis/train/testResults" << called << ".txt";
    string filename = fnstream.str();
    ofstream trainingOutput(filename.c_str(),ios::out);
    cout << "trainingSetError: " << details(trainingData, trainingOutput);
    trainingOutput.close();

    fnstream.str("");
    fnstream << "analysis/test/testResults" << called << ".txt";
    filename = fnstream.str();
    ofstream testOutput(filename.c_str(),ios::out);
    cout << " testSetError: " << details(testData, testOutput);
    testOutput.close();
    cout << " frame #" << called << endl;
    
    called++;
}

double MultiLayerPerceptron::details(SurfData* data, std::ostream& os)
{
    if (!data) {
	    cerr << "No data in details" << endl;
    }
    const SurfPoint* currentPoint;
    double output;
    double actual;
    double error = 0.0;
    for (int i = 0; i < data->size(); i++) {
        currentPoint = data->getPoint(i);
 	actual = currentPoint->getF(responseIndex);	
	output = onePattern(currentPoint->getX());
	currentPoint->write(os);
	os << setw(15) << output
		<< setw(15) << actual - output
		<< endl;
	error += (actual - output) * (actual - output);
    }
    return .5 * error;
	
}
    

bool MultiLayerPerceptron::makeFrameThisIteration(int iteration)
{
    return false;
    return (
	     (iteration < 100 && iteration % 10 == 0) || (iteration % 20 == 0)
	     //(iteration < 10000 && iteration % 25 == 0) ||
	     //(iteration < 10000 && iteration % 1000 ==0) 
	   );
}

void MultiLayerPerceptron::setTestData(SurfData* testData)
{
    this->testData = testData;
}
//double MultiLayerPerceptron::baselineError()
//{
//    double sum = 0;
//    for (int i = 0; i < numPatterns; i++) {
//	sum += patterns[i].target;
//    }
//    double avg = sum / static_cast<double>(numPatterns);
//    double sse = 0;
//    for (int j = 0; j < numPatterns; j++) {
//	sse += pow(patterns[j].target - avg, 2);
//    }
//    return .5 * sse;    
//
//}
	
//double MultiLayerPerceptron::rosenbrock(double x1, double x2)
//{
//    //return x1*x1 + x2*x2;
//    return x1*x1 + x2*x2 - 10.0*cos(2.0*PI*x1) - 10.0*cos(2.0*PI*x2) + 20.0;
//    //return (100.0*(x1 - x2*x2)*(x1-x2*x2) + (1-x2)*(1-x2)) ;
//}


//void MultiLayerPerceptron::interpolate()
//{
//    string filename;
//    ostringstream makestring;
//    ofstream totalErrors("Interpolations2/totalErrors.txt",ios::out);
//    double interpolationError = 0.0;
//    for (int i = 0; i < numPatterns - 1; i++) {
//	for (int j = i + 1; j < numPatterns; j++) { 
//	    interpolationError = 0.0;
//	    makestring << "Interpolations2/interp" << i << "_" << j << ".txt";
//	    filename = makestring.str();
//	    cout << filename << endl;
//	    makestring.str("");
//	    ofstream outfile(filename.c_str(),ios::out);
//	    double r = 0.0;
//	    double step = 0.01;
//	    double squaredError = 0.0;
//	    Pattern newp;
//	    Pattern& first = patterns[i];
//	    Pattern& second = patterns[j];
//	    while (r < 1.01) {
//		newp.inputs[0] = (1.0-r)*first.inputs[0] +  r*second.inputs[0];
//		newp.inputs[1] = (1.0-r)*first.inputs[1] +  r*second.inputs[1];
//	        newp.target = rosenbrock(newp.inputs[0], newp.inputs[1]);
//		pattern = &newp;
//		fire();
//		squaredError = .5 * pow( outputs2[0] - pattern->target, 2);
//		interpolationError += squaredError;
//	 	outfile << setprecision(10)
//	                << setw(20) << newp.inputs[0]
//			<< setw(20) << newp.inputs[1]
//		        << setw(20) << r
//			<< setw(20) << newp.target
//			<< setw(20) << outputs2[0]
//			<< setw(20) << squaredError
//			<< endl;
//		r += step;
//	    }
//	    outfile.close();
//	    totalErrors << setw(20) << interpolationError << "\t"
//			<< setw(40) << filename << endl;
//	}
//    }
//    totalErrors.close();
//
//
//}

