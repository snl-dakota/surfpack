// These lines were automatically prepended
#ifndef SURFPACK_CONFIG_H
#define SURFPACK_CONFIG_H
#include "config.h"
#endif
// End of prepended lines

#include "unittests.h"

using namespace std;

const string fullPath(const string filename)
{
  return dpath + "/" + filename;
}
