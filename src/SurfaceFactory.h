// These lines were automatically prepended
#ifndef SURFPACK_CONFIG_H
#define SURFPACK_CONFIG_H
#include "config.h"
#endif
// End of prepended lines

#ifndef __SURFPACK_FACTORY_H__
#define __SURFPACK_FACTORY_H__

class Surface;
class SurfData;
#include <vector>
#include "SurfpackParser.h"


namespace SurfaceFactory {
  Surface* createSurface(const std::string& filename);
  Surface* createSurface(const std::string& type, SurfData* surfData);
  //Surface* createSurface(const std::string*, SurfData* sd, 
  //  const SurfpackParser::ArgList& arglist);
  //Surface* createSurface(const std::string& type, SurfData* surfData, unsigned order);
}

#endif
