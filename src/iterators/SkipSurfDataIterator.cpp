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

#include "SkipSurfDataIterator.h"

using namespace std;

// constructors/destructor
////////////////////////////////////

// copy constructor
//
//SkipSurfDataIterator::SkipSurfDataIterator(const SkipSurfDataIterator & iterator) 
//{ 
//} 
 

/// iterates over surfData, skipping specified points
SkipSurfDataIterator::SkipSurfDataIterator(SurfData* surfData) : AbstractSurfDataIterator(surfData)
{
}

SkipSurfDataIterator::~SkipSurfDataIterator() 
{

}

// public methods
////////////////////////////////////

/// get the next element in the iterator
SurfPoint * SkipSurfDataIterator::nextElement() 
{ 
   if(!surfData) {
	  cerr << "Null SurfData in SkipSurfDataIterator::nextElement()" << endl;
	  return 0;
   } else if (currentIndex < surfData->size()) {
       ++currentIndex; 
   }
   return currentElement(); 
}

/// get the current element from the iterator
SurfPoint * SkipSurfDataIterator::currentElement() 
{ 
   if(!surfData) {
	  cerr << "Null SurfData in SkipSurfDataIterator::currentElement()" << endl;
	  return 0;
   }
   while(currentIndex < surfData->size() && skipIndices.find(order[currentIndex]) !=skipIndices.end() ) {
	   ++currentIndex;
   }
   if (currentIndex >= surfData->size()) {
       return 0; 
   } else {
       return surfData->getPoint(order[currentIndex]);
   }
}

/// is the iterator at the end
bool SkipSurfDataIterator::isEnd() 
{ 
    currentElement();
    return currentIndex >= surfData->size();
}


/// skip over point at <index> during iteration
void SkipSurfDataIterator::skipPoint(int index)
{ 
   if (index >= 0 && index < surfData->size()) {
       skipIndices.insert(index); 
   }
}

/// return the number of elements in one complete iteration 
int SkipSurfDataIterator::getElementCount() const
{
    return surfData->size() - skipIndices.size();
}

/// unmark all elements that are currently skipped
void SkipSurfDataIterator::unSkipAll()
{
    skipIndices.clear();
}
