// ----------------------------------------------------------
// Project: SURFPACK++
//
// File:        SurfSurfDataIterator.h
// Author:      Mark Richards 
// Created:     August 12, 2003
// Modified:
//
// Description: 
// + SurfSurfSurfDataIterator class - Iterator for the SurfData class
// ----------------------------------------------------------

#ifndef __DATA_ITERATOR_H__
#define __DATA_ITERATOR_H__

#include <vector>

#include "AbstractSurfDataIterator.h"

class SurfDataIterator : public AbstractSurfDataIterator
{
private:

// private data members
////////////////////////////////////
   


// private methods
////////////////////////////////////

   // copy constructor
   //
//   SurfDataIterator(const SurfDataIterator &); 

   // assignment operator
   //
//   SurfDataIterator & operator=(const SurfDataIterator &);

   SurfDataIterator() {}


public:

// constructors/destructor
////////////////////////////////////
 
   /// default constructor
   SurfDataIterator(SurfData* surfData);

   virtual ~SurfDataIterator();

// public methods
////////////////////////////////////

   /// return the next element in the iterator
   virtual SurfPoint * nextElement();

   /// return the current element from the iterator
   virtual SurfPoint * currentElement();

   /// get the number of elements in this iterator 
   virtual int getElementCount() const;
};

#endif
