// These lines were automatically prepended
#ifndef SURFPACK_CONFIG_H
#define SURFPACK_CONFIG_H
#include "config.h"
#endif
// End of prepended lines

#ifndef UNITTESTS_H
#define UNITTESTS_H

#include <string>


// DATADIR is defined in configure.ac as $(top_srcdir)/share
const std::string dpath = DATADIR;
const std::string fullPath(const std::string filename);
const unsigned unsignedZero = 0;
const double doubleZero = 0.0;
const int intZero = 0;


#endif
