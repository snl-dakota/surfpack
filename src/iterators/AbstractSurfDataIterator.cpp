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

#include <algorithm>
#include <vector>
#include <iostream>
#include "SurfData.h"
#include "AbstractSurfDataIterator.h"

using namespace std;

// constructors/destructor
////////////////////////////////////
 
/// iterate over elements in surfData
AbstractSurfDataIterator::AbstractSurfDataIterator(SurfData& sd, unsigned response_index)  
  : sd(sd), response_index(response_index)
{
    sort();
}


// copy constructor
//AbstractSurfDataIterator::AbstractSurfDataIterator(const AbstractSurfDataIterator & di) 
//{ 
//    this->sd = di.sd;
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
    return currentIndex >= sd.size(); 
}

/// set the iterator to the first element 
void AbstractSurfDataIterator::toFront() 
{ 
    currentIndex = 0;
}

void AbstractSurfDataIterator::sort()
{
    order.resize(sd.size());
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

unsigned AbstractSurfDataIterator::xSize() const
{
    return sd.xSize();
}

/// Return index of response variable associated with this iterator
unsigned AbstractSurfDataIterator::responseIndex() const
{
  return response_index;
}
