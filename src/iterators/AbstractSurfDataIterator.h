// ----------------------------------------------------------
// Project: SURFPACK++
//
// File:        AbstractSurfDataIterator.h
// Author:      Mark Richards 
// Created:     August 12, 2003
//
// Description: 
// + AbstractSurfDataIterator class - Iterator for the AbstractData class
// ----------------------------------------------------------

#ifndef __ABSTRACT_SURFDATA_ITERATOR_H__
#define __ABSTRACT_SURFDATA_ITERATOR_H__

class SurfPoint;
#include "SurfData.h"
class AbstractSurfDataIterator
{
private:

// private data methods
////////////////////////////////////

   // copy constructor
   //
//   AbstractSurfDataIterator(const AbstractSurfDataIterator &); 

   // assignment operator
   //
//   AbstractSurfDataIterator & operator=(const AbstractSurfDataIterator &);

protected:
   AbstractSurfDataIterator(); 

// data members
////////////////////////////////////

   SurfData& sd;
   unsigned response_index;
   unsigned currentIndex;
   std::vector<unsigned> order;

public:

// constructors/destructor
////////////////////////////////////
 
   /// The SurfData object which will be iterated over must be provided in the constructor
   AbstractSurfDataIterator(SurfData& sd, unsigned response_index = 0);
   
   virtual ~AbstractSurfDataIterator(); 

// public methods
////////////////////////////////////

   /// return the next element in the iterator
   virtual SurfPoint& nextElement() = 0;

   /// return the current element from the iterator
   virtual SurfPoint& currentElement() = 0;

   /// reset the iterator to the first element 
   virtual void toFront();

   /// is the iterator at the end
   virtual bool isEnd();

   /// return the number of SurfPoints that comprise a full iteration
   virtual unsigned elementCount() const = 0;

   /// return the number of dimensions in the SurfData
   unsigned xSize() const;

   /// Return index of response variable associated with this iterator
   unsigned responseIndex() const;

   /// produce a randomly ordered iteration sequence
   virtual void shuffle();

   /// restore sequence of SurfPoints to its original order
   virtual void sort(); 
};

#endif
