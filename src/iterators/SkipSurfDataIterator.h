// ----------------------------------------------------------
// Project: SURFPACK++
//
// File:        SkipSurfDataIterator.h
// Author:      Mark Richards 
// Created:     August 12, 2003
// Modified:
//
// Description: 
// + SkipSurfDataIterator class - Iterator for the SurfData class
// ----------------------------------------------------------


#ifndef __SKIP_SURFDATA_ITERATOR_H__
#define __SKIP_SURFDATA_ITERATOR_H__

#include <set>

#include "AbstractSurfDataIterator.h"

class SkipSurfDataIterator : public AbstractSurfDataIterator
{
private:

// private data members
////////////////////////////////////

   std::set<int> skipIndices;

// private methods
////////////////////////////////////

//   SkipSurfDataIterator(const SkipSurfDataIterator &); 
//   SkipSurfDataIterator & operator=(const SkipSurfDataIterator &);

   // default constructor
   //
   SkipSurfDataIterator() {}

public:

// constructors/destructor
////////////////////////////////////
 

   SkipSurfDataIterator(SurfData& sd, unsigned response_index = 0);

   virtual ~SkipSurfDataIterator();

// public methods
////////////////////////////////////

   /// get the next element in the iterator
   virtual SurfPoint& nextElement();

   /// get the current element from the iterator
   virtual SurfPoint& currentElement();

   /// is the iterator at the end 
   virtual bool isEnd();

   /// return number of elements in one complete iteration 
   virtual unsigned elementCount() const;
   
   /// add one point only to the list of points to skip
   void skipPoint(unsigned index);

   /// unmark all elements currently being skipped
   void unSkipAll();
 
};

#endif
