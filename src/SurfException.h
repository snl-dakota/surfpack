// ----------------------------------------------------------
// Project: SURFPACK++
//
// File:        SurfException.h
// Author:      Eric Cyr
// Created:     June 7, 2002
// Modified:    June 24, 2002
//
// Description: 
// + SurfException class - A very generic exception class,
//   it only contains a simple string message and 
//   an integer error code.
// ----------------------------------------------------------

#ifndef __SURF_EXCEPTION_H__
#define __SURF_EXCEPTION_H__

#include <iostream>
#include <string>

class SurfException
{
private:

// Private data memebers
////////////////////////////////////

  std::string message;  // the message for this exception
  int code;        // the error code for this exception

public:

   // Constructs a SurfException object with some meaningless
   // defaults
   //  
   SurfException() 
      : message(""), code(0) {}

   // Constructs a SurfException object with a set message
   // 
   SurfException(const std::string & m)
      : message(m), code(0) {}

   // Constructs a SurfException object with a set message
   // and error code 
   //
   SurfException(const std::string & m,int c)
      : message(m), code(c) {}

   // do nothing destructor
   //
   virtual ~SurfException() {}

   // Get the message associated with this exception
   //
   std::string getMessage() { return message; }

   // Get the error code associated with this exception
   //
   int getErrorCode() { return code; }

   // mechanism to help printing
   //
   std::ostream & print(std::ostream & os) const
   { return os << message << " " << code; }
};

// implmentation for this is in SurfData.cpp
//
std::ostream & operator<<(std::ostream &,const SurfException &);


#endif
