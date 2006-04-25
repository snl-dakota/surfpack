/*  _______________________________________________________________________

    DAKOTA: Design Analysis Kit for Optimization and Terascale Applications
    Copyright (c) 2001, Sandia National Laboratories.

    Surfpack: A Software Library of Multidimensional Surface Fitting Methods

    Surfpack is distributed under the DAKOTA GNU General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */
#ifndef FLEX_WRAPPER_H
#define FLEX_WRAPPER_H
#ifdef HAVE_CONFIG_H
#include "surfpack_config.h"
#endif
#include "surfpack_system_headers.h"
class FlexWrapper
{
public:
  FlexWrapper();
  ~FlexWrapper();
  void setParseStreams(const std::string* input_string, const std::string* output_string);
  int nextToken();
  const char* currentToken();
private:
  FILE* infile;
  FILE* outfile;
};

#endif
