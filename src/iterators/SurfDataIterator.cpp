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

#include <vector>
#include <algorithm>
#include <iostream>
#include "SurfDataIterator.h"
#include "SurfPoint.h"
#include "SurfData.h"

using namespace std;

// constructors/destructor
////////////////////////////////////
 
/// iterate over elements in surfData
SurfDataIterator::SurfDataIterator(SurfData& sd, unsigned response_index) 
  : AbstractSurfDataIterator(sd, response_index)
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
SurfPoint& SurfDataIterator::nextElement() 
{ 
  if (currentIndex < sd.size()) {
      currentIndex++;
  }
  return currentElement();
}

/// return the current element from the iterator
SurfPoint& SurfDataIterator::currentElement() 
{ 
  return sd.Point(order[currentIndex]);
}

/// return the number of elements in one complete iteration
unsigned SurfDataIterator::elementCount() const 
{ 
  return sd.size();
}    
