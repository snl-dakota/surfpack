// ----------------------------------------------------------
// Project: SURFPACK++
//
// File:        AbstractSurfDataIterator.cpp
// Author:      Mark Richards 
// Created:     August 12, 2003
// Modified:
//
// Description: 
// + SurfAbstractSurfDataIterator class - Iterator for the SurfData class
// ----------------------------------------------------------

#include "AbstractSurfDataIterator.h"

using namespace std;

// constructors/destructor
////////////////////////////////////
 
/// iterate over elements in surfData
AbstractSurfDataIterator::AbstractSurfDataIterator(SurfData* surfData)  : surfData(surfData)
{
    sort();
}


// copy constructor
//AbstractSurfDataIterator::AbstractSurfDataIterator(const AbstractSurfDataIterator & di) 
//{ 
//    this->surfData = di.surfData;
//    this->order = di.order;
//    this->currentIndex = di.currentIndex;
//}

// destructor
AbstractSurfDataIterator::~AbstractSurfDataIterator() 
{

}

// public methods
////////////////////////////////////


/// is the iterator at the end
bool AbstractSurfDataIterator::isEnd() 
{ 
    return currentIndex >= surfData->size(); 
}

/// set the iterator to the first element 
void AbstractSurfDataIterator::toFront() 
{ 
    currentIndex = 0;
}

void AbstractSurfDataIterator::sort()
{
    order.resize(surfData->size());
    for (unsigned i = 0; i < order.size(); i++) {
	order[i] = i;
    }
    currentIndex = 0;
}

void AbstractSurfDataIterator::shuffle()
{
    random_shuffle(order.begin(),order.end());
    currentIndex = 0;
}

int AbstractSurfDataIterator::getDimension() const
{
    return surfData->getDimension();
}

