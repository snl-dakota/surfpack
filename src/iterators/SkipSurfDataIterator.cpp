// Project: SURFPACK++
//
// File:        SkipSurfDataIterator.cpp
// Author:      Mark Richards 
// Created:     August 12, 2003
// Modified:
//
// Description: 
// + SkipSurfDataIterator class - Iterator for the SurfData class
// ----------------------------------------------------------

#include <vector>
#include <algorithm>
#include <iostream>
#include "SkipSurfDataIterator.h"
#include "SurfPoint.h"
#include "SurfData.h"

using namespace std;

// constructors/destructor
////////////////////////////////////

// copy constructor
//
//SkipSurfDataIterator::SkipSurfDataIterator(const SkipSurfDataIterator & iterator) 
//{ 
//} 
 

/// iterates over surfData, skipping specified points
SkipSurfDataIterator::SkipSurfDataIterator(SurfData& sd, unsigned response_index) 
  : AbstractSurfDataIterator(sd, response_index)
{
}

SkipSurfDataIterator::~SkipSurfDataIterator() 
{

}

// public methods
////////////////////////////////////

/// get the next element in the iterator
SurfPoint& SkipSurfDataIterator::nextElement() 
{ 
  if (currentIndex < sd.size()) {
      ++currentIndex; 
  }
  return currentElement(); 
}

/// get the current element from the iterator
SurfPoint& SkipSurfDataIterator::currentElement() 
{ 
  while(currentIndex < sd.size() && skipIndices.find(order[currentIndex]) !=skipIndices.end() ) {
          ++currentIndex;
  }
  return sd.Point(order[currentIndex]);
   
}

/// is the iterator at the end
bool SkipSurfDataIterator::isEnd() 
{ 
  currentElement();
  return currentIndex >= sd.size();
}


/// skip over point at <index> during iteration
void SkipSurfDataIterator::skipPoint(unsigned index)
{ 
  if (index >= 0 && index < sd.size()) {
    skipIndices.insert(index); 
  }
}

/// return the number of elements in one complete iteration 
unsigned SkipSurfDataIterator::elementCount() const
{
  return sd.size() - skipIndices.size();
}

/// unmark all elements that are currently skipped
void SkipSurfDataIterator::unSkipAll()
{
  skipIndices.clear();
}
