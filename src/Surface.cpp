// ----------------------------------------------------------
// Project: SURFPACK++
//
// File:        Surface.cpp
// Author:      Tony Giunta
// Created:     June 13, 2002
// Modified:    Eric Cyr - June 18, 2002        
//                         June 20, 2002
//                         July 15-16, 2002        
//
// Description: 
// + The Surface class provides an interface for the
//   different surface types, this class cannot be instantiated
// ----------------------------------------------------------

#include "Surface.h"
#include "SurfData.h"
#include <fstream>
#include <iostream>
#include <iomanip>

using namespace std;

#ifdef __TESTING_MODE__
   int Surface::constructCount;
   int Surface::destructCount;
#endif

// destructor, just removes itself from any surfaces
//
Surface::~Surface() 
{ 
   if(surfData)
      surfData->removeListener(this);

//   delete surfDataIter;

#ifdef __TESTING_MODE__
   destructCount++;
#endif
}

// a method to print the surface info to
//
// os: the stream to write to
// return: returns the same stream passed in
//
//ostream & Surface::print(ostream & os = cout) 
//{
//   os << "Surface type = " << surfType << endl;
//   for(unsigned int i=0; i<coefficients.size(); i++ ) {
//	 os << coefficients[i] << endl;
//   }
//   return os;
//}

//ostream & Surface::printClean(ostream & os = cout) throw(SurfException)
//{
//   static double tolerance = 1.0e-10;
//   os << "Surface type = " << surfType << endl;
//   for(unsigned int i=0; i<coefficients.size(); i++ ) {
// // really small numbers are effectively zero; the printout is more visually accessible if
// // these numbers show up as zeroes.
//      if (coefficients[i] < tolerance && -coefficients[i] < tolerance) {
//	 os << 0 << endl;
//      } else {
//	 os << coefficients[i] << endl;
//      }
//   }
//   return os;
//}

// basic constructor
//
// s: type of surface being created
//
//Surface::Surface(const string & s) :
//   surfType(s), surfData(0), needsRebuild(true), responseIndex(-1) 
//{
//#ifdef __TESTING_MODE__
//   constructCount++;
//#endif
//
//   surfDataIter = new DataIterator();
//}

// constructor that takes the coefficients for this object
//
// s: type of surface being created
// coeff: coefficients for this surface 
//
//Surface::Surface(const vector<double> & coeff) :
//   surfData(0), needsRebuild(true), responseIndex(-1)
//{
//#ifdef __TESTING_MODE__
//   constructCount++;
//#endif
//
////   surfDataIter = new DataIterator();
//}

// constructor that takes the data set for this object
// adds itself as a listener to the SurfData object
//
// s: type of surface being created
// data: the SurfData object associated with this surface 
//
Surface::Surface(SurfData * surfData, int responseIndex) :
   surfData(surfData), needsRebuild(surfData != 0), responseIndex(responseIndex) 
{ 
#ifdef __TESTING_MODE__
   constructCount++;
#endif

   if(surfData) 
      surfData->addListener(this); 

   //surfDataIter = new DataIterator();
}

// constructor that takes the coefficients for this object
//
// s: type of surface being created
// coeff: coefficients for this surface 
//
//Surface::Surface(const string & s,SurfData * data,const vector<double> & coeff) :
//   surfType(s), surfData(data), coefficients(coeff), needsRebuild(true), responseIndex(-1)
//{ 
//#ifdef __TESTING_MODE__
//   constructCount++;
//#endif
//
//   if(surfData)
//      surfData->addListener(this);  
//
//   surfDataIter = new DataIterator();
//}

// basic copy constructor, ignores most things and copies only what
// it needs
//
// surf: the surface to be copied
//
//Surface::Surface(const Surface & surf) : 
//   surfType(surf.surfType), surfData(0), 
//   coefficients(surf.coefficients), needsRebuild(true), responseIndex(-1)
//{ 
//#ifdef __TESTING_MODE__
//   constructCount++;
//#endif
//
//   surfDataIter = new DataIterator();
//}


ostream & operator<<(ostream & os,Surface & surface)
{ return surface.write(os); }

// set the index of the surface to be modeled
//
// index: index of the surface, make sure the index
//        is within range
//
bool Surface::setResponseIndex(int index)
{ 
   if(surfData==0)
   {
      //responseIndex = index;
      return false;
   }

   int respCnt = surfData->getResponseCount();

   if (index < respCnt) {
       responseIndex = index;
       return true;
   } else {
       cerr << "Response index " << index << " is not valid for this data" << endl;
   }
   // here I'm allowing the user to have a data set with
   // no responses but in anticipation I allow them to
   // use the first spot
   //if(index>=0 && (index<respCnt || (respCnt<=0 && index==0))) 
   //{
   //   responseIndex = index; 
   //   needsRebuild = true;
 
   //   return true;
   //}

   // index out of range 
   return false;
}

int Surface::getResponseIndex()
{ return responseIndex; }

// set the iterator for this class 
// uses the clone funtion to copy the object
// copy constructors are not polymorphic
//
//void Surface::setIterator(const AbstractDataIterator & iter)
//{
//   AbstractDataIterator * delMe = surfDataIter;
//   surfDataIter = iter.clone();
//   delete delMe; 
//   needsRebuild = true;
//}

//void Surface::setData(SurfData * data)
//{
//   // clean up all the old sutff
//   if(surfData)
//      surfData->removeListener(this);
//
//   // set to the new one
//   surfData = data;
//  
//   // setup the new one
//   if(surfData)
//      surfData->addListener(this);
//}

double Surface::evaluate(const std::vector<double> & x) throw(SurfException)
{
   // check to be sure everything is OK
   if(needsRebuild) {
       build();
   }

   //if (x.size() != static_cast<unsigned>(getDimension())) {
   //    cerr << "Cannot evaluate point: wrong number of dimensions" << endl;
   //    return 0.0;
   //}
   
   return calculate(x);
}

void Surface::evaluate(SurfData* surfData)
{
    int newindex = surfData->addResponse();
    SurfDataIterator surfIt(surfData);
    surfIt.toFront();
    while (!surfIt.isEnd()) {
	SurfPoint* current = surfIt.currentElement();
        double estimatedResponse = evaluate(current->getX());
	current->setF(newindex, estimatedResponse);
	surfIt.nextElement();
    }
}

/// Evaluate the empirical model at the points in surfData and output
/// the points and their evaluations to os
double Surface::test(SurfData* surfData, ostream & os)
{
    if (!surfData) {
	    cerr << "Null SurfData object in evaluate" << endl;
	    return 0.0;
    }
    ofstream outfile("results.txt", ios::out);
    double error = 0.0;
    for (int i = 0; i < surfData->size(); i++) {
	SurfPoint* currentPoint = surfData->getPoint(i);
	if (!currentPoint) {
		cerr << "Null pointer returned in Surface::evaluate(SurfData)" << endl;
	}
	double estimatedResponse = evaluate(currentPoint->getX());

	// print out the points and deviations to a file
        currentPoint->write(outfile);
	outfile << setw(20) << estimatedResponse 
	        << setw(20) << (currentPoint->getF() - estimatedResponse)
		<< endl;	
	error += (currentPoint->getF() - estimatedResponse) * (currentPoint->getF() - estimatedResponse);
    }
    error *= .5;
    outfile.close();
    return error;
}    


//double Surface::evaluate(double * x,int dim) throw(SurfException)
//{
//   vector<double> xVec;
//
//   // convert the array to an STL vector
//   for(int i=0;i<dim;i++)
//      xVec.push_back(x[i]);
//
//   return evaluate(xVec);
//}

//double Surface::evaluate(const vector<double> & x,vector<double> & grad)
//   throw(SurfException)
//{
//   // check to be sure everything is OK
//   if(surfData)
//   {
//      surfData->getIterator(surfDataIter);
//
//      // make sure the domainis the right size only if you have a
//      // SurfData object associated with you
//      if((int) x.size() < getDimension()) 
//         throw SurfException("Point to evaluate at is the wrong dimension");
//
//      if(surfDataIter->getElementCount() < getMinPointCount()) 
//         throw SurfException("Point to evaluate at is the wrong dimension");
//   }
//
//   if(needsRebuild) 
//      build();
//
//   return calculate(x,grad);
//}

//double Surface::evaluate(const double * x,double * grad,int dim) throw(SurfException)
//{
//   vector<double> xVec(dim);
//   vector<double> gradVec(dim);
//
//   // convert into an stl vector
//   for(int i=0;i<dim;i++)
//      xVec[i] = x[i];
//
//   double answer = evaluate(xVec,gradVec);
//
//   if((int) gradVec.size()!=dim) 
//      throw SurfException("Gradient is the incorrect size");
//
//   // convert back into a array
//   for(int i=0;i<dim;i++)
//      grad[i] = gradVec[i];
//
//   return answer;
//}

//void Surface::evaluateGrad(const std::vector<double> & x,std::vector<double> & grad)
//   throw(SurfException)
//{
//   // check to be sure everything is OK
//   if(surfData)
//   {
//      surfData->getIterator(surfDataIter);
//
//      // make sure the domainis the right size only if you have a
//      // SurfData object associated with you
//      if((int) x.size() < getDimension()) 
//         throw SurfException("Point to evaluate at is the wrong dimension");
//
//      if(surfDataIter->getElementCount() < getMinPointCount()) 
//         throw SurfException("Point to evaluate at is the wrong dimension");
//   }
//
//   if(needsRebuild) 
//      build();
//
//   calculate(x,grad);
//}
//
//void Surface::evaluateGrad(const double * x,double * grad,int xGradSz)
//   throw(SurfException)
//{
//   vector<double> gradVec;
//   vector<double> xVec;
//
//   // convert the array to STL vector
//   for(int i=0;i<xGradSz;i++)
//      xVec.push_back(x[i]);
//
//   evaluateGrad(xVec,gradVec);
//
//   // make sure everything is the right size (this singals in
//   // error internally)
//   if((int) gradVec.size()!=xGradSz) 
//      throw SurfException("Gradient is the incorrect size");
//
//   for(unsigned int i=0; i<(unsigned int) xGradSz && i<gradVec.size();i++)
//      grad[i] = gradVec[i];
//}

void Surface::build() throw(SurfException)
{
    SurfDataIterator iter(surfData);
    build(&iter);
}

void Surface::build(AbstractSurfDataIterator* iter) throw(SurfException)
{
   if(surfData) {
       if (iter->getElementCount() < getMinPointCount()) {
           cerr << "Not enough points to create this surface" << endl;
       } else {
           calculateInternals(iter);
           needsRebuild = false;
       }
   } else {
       cerr << "No SurfData: surface cannot be built" << endl;
   }
}

//double Surface::calculate(const vector<double> & x) throw(SurfException)
//{
//   vector<double> falseGrad;
//   falseGrad.resize(x.size());
//   double result = calculate(x,falseGrad,false);
//   falseGrad.resize(4);
//   return result; 
//}

double Surface::computePressStatistic()
{
   double fpress = 0;
   //DataIterator dataIterator;
   //surfData->getIterator(dataIterator);
   const SurfPoint* currentPoint = 0;
   int i = 0;
   SkipSurfDataIterator skipIt(surfData);
   while (i < surfData->size()) {
	
	// leave out one point from the data set
	currentPoint = surfData->getPoint(i); 
	assert(currentPoint);
	skipIt.unSkipAll();
	skipIt.skipPoint(i);
	//cout << "Skipping: " << i << endl;
	build(&skipIt);
	double fTilde = evaluate(currentPoint->getX());
	
	// actual value of the function at this point
	double fx = currentPoint->getF();

	// add the square of the deviation fTilde - fx
	//cout << "Predicted: " << fTilde << " Observed: " << fx << " Difference: " << fx - fTilde << endl;
	fpress += (fTilde - fx) * (fTilde - fx);
	//dataIterator.nextElement();
	i++;
	
   }
   skipIt.unSkipAll();
   build(&skipIt);
   fpress = sqrt(fpress/i);
   return fpress;
}
   	     	
//double Surface::getPressStatistic()
//{
//	if (needsRebuild) {
//		build();
//		return computePressStatistic();
//	}
//
//	return pressStatistic;
//}

bool significant(double val)
{
	static double tolerance = 1.0e-10;
	return val > tolerance || -val > tolerance;
}
	
/// Save the surface to a file in binary format
void Surface::saveBinary(std::string filename)
{

}

/// Load the surface from a file
void Surface::loadBinary(std::string filename)
{

}

double Surface::errorMetric(string metricName)
{
	if (metricName == "press") {
		cout << "Call Press from Surface" << endl;
		return computePressStatistic();
	//} else if (metricName == "rsquared") {
	//	cout << "Call rsquared from Surface" << endl;
	//	return rsquared();
	} else {
		cout << "No error metric of that type in this class" << endl;
	}
	return 0;

}

double Surface::press()
{
	cout << "Surface::press" << endl;
	return 0;

}

double Surface::rsquared()
{
	cout << "Surface::rsquared" << endl;
	return 0;

}
