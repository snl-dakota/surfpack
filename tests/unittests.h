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

void writePoint1Files();
void writePoint2Files();
void writeRastriginAndClaimsTooManyFiles();
void writeManyPtsFiles();
void writeOneDimQuadratic();
void writeUnknownSurfaceFile();
void setOstreamFlags(std::ostream& os);
void initialize();
void cleanup();
bool matches(double observed, double target, double margin = 1.0e-2);

#endif
