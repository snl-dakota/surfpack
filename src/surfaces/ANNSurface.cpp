#include "ParticleSwarmMLP.h"
#include "BackpropMLP.h"
#include "ANNSurface.h"
#include "SurfData.h"
#include <iostream>
#include <cmath>

using namespace std;

//
//ANNSurface::ANNSurface(const std::vector<double> & coeff)
//   : Surface("ANNSurface",coeff) { needsRebuild = false; }
//
//ANNSurface::ANNSurface(SurfData * data)
//   : Surface("ANNSurface",data) { }
//
//ANNSurface::ANNSurface(const ANNSurface & surf)
//   : Surface(surf) { needsRebuild = false; }
//
ANNSurface::~ANNSurface() 
{
    if (ownMLP) {
        delete mlp;
    }
}

ANNSurface::ANNSurface(SurfData* surfData, int responseIndex, int size, int maxTrainingIterations, double goalTrainingError, int trainingMethodCode)
    : Surface(surfData, responseIndex) 
{
    if (trainingMethodCode == MultiLayerPerceptron::particleSwarm) {
	    mlp = new ParticleSwarmMLP(surfData,responseIndex,size,maxTrainingIterations,goalTrainingError);
    } else {
	    mlp = new BackpropMLP(surfData,responseIndex,size,maxTrainingIterations,goalTrainingError);
    }
    ownMLP = true;
}

ANNSurface::ANNSurface(SurfData* surfData, MultiLayerPerceptron* mlp, int responseIndex)
    : Surface(surfData, responseIndex), mlp(mlp)
{
    ownMLP = false;
    if (!mlp) {
	cerr << "Null mlp in ANNSurface constructor" << endl;
    }
}

int ANNSurface::getMinPointCount(int val) const
{ 
    return val*2;
}

int ANNSurface::getDimension() const
{
    return mlp->getNumInputs();
}
   
double ANNSurface::calculate(const std::vector<double> & x) throw(SurfException)
{ 
    return mlp->onePattern(x);
}
    
void ANNSurface::calculateInternals(AbstractSurfDataIterator* iter)
{
    mlp->trainNetwork();
}

ostream & ANNSurface::write(ostream & os) 
{
    os << "ANNSurface cannot currently be written" << endl;
    return os;
}

string ANNSurface::getType() const
{
    return "Artificial Neural Network";
}

void ANNSurface::save(string filename)
{
    mlp->save(filename);
}
