// ----------------------------------------------------------
// Project: SURFPACK++
//
// File:        ANNSurface.h
// Author:      Mark Richards	 
// Created:     July 18, 2003
// Modified:    
//
// Description: 
// + The ANNSurface class uses an artificial neural network to fit the points 
// ----------------------------------------------------------

#ifndef __ANN_SURFACE_H__
#define __ANN_SURFACE_H__

#include <iostream>
#include <vector>
#include <string>

#include "SurfException.h"
#include "SurfPoint.h"
#include "Surface.h"
#include "SurfDataIterator.h"
#include "SkipSurfDataIterator.h"


class MultiLayerPerceptron;

class ANNSurface : public Surface
{
private:
   ANNSurface() {}

public:

// constructor/destructor
////////////////////////////////

//   ANNSurface(const std::vector<double> &);
//   ANNSurface(SurfData *);
//   ANNSurface(SurfData *,const std::vector<double> &);
//   ANNSurface(const ANNSurface &);
     ANNSurface(SurfData* surfData, int responseIndex = 0, int size = 10, int maxTrainingIterations = 10000, double goalTrainingError = 1e-8, int trainingMethodCode = 1); 
     ANNSurface(SurfData* surfData, MultiLayerPerceptron* mlp, int responseIndex = 0);
     ANNSurface(std::istream& is);
   
   // a do nothing destructor
   //
   virtual ~ANNSurface(); 

// member functions
////////////////////////////////

   virtual int getMinPointCount(int dim) const;
   virtual int getDimension() const;
   virtual std::ostream & write(ostream & os = cout); 
   virtual void save(std::string filename);
   virtual std::string getType() const;
   

protected:
   MultiLayerPerceptron* mlp;
   bool ownMLP;

   virtual double calculate(const std::vector<double> & x) throw(SurfException);
   virtual void calculateInternals(AbstractSurfDataIterator* iter);
};

#endif
