// ----------------------------------------------------------
// Project: SURFPACK++
//
// File:        SurfDataIterator.cpp
// Author:      Mark Richards 
// Created:     August 12, 2003
// Modified:
//
// Description: 
// + SurfSurfDataIterator class - Iterator for the SurfData class
// ----------------------------------------------------------

#include "SurfDataIterator.h"

using namespace std;

// constructors/destructor
////////////////////////////////////
 
/// iterate over elements in surfData
SurfDataIterator::SurfDataIterator(SurfData* surfData) : AbstractSurfDataIterator(surfData)
{
}

// copy constructor
//SurfDataIterator::SurfDataIterator(const SurfDataIterator & di) 
//{ 
//}

// destructor
// 
SurfDataIterator::~SurfDataIterator() { }

// public methods
////////////////////////////////////

/// return the next element in the iterator
SurfPoint * SurfDataIterator::nextElement() 
{ 
    if (currentIndex < surfData->size()) {
	currentIndex++;
    }
    return currentElement();
}

/// return the current element from the iterator
SurfPoint * SurfDataIterator::currentElement() 
{ 
    if (currentIndex < surfData->size()) { 
	return surfData->getPoint(order[currentIndex]);
    } 
    return 0; 
}

/// return the number of elements in one complete iteration
int SurfDataIterator::getElementCount() const 
{ 
    return surfData->size();
}    
