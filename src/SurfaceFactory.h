#include "config.h"

#ifndef __SURFPACK_FACTORY_H__
#define __SURFPACK_FACTORY_H__

#include <vector>

#include "SurfpackParser.h"

class Surface;
class SurfData;
/// The createSurface methods are intended to be a sort of virtual constructor.
/// When new Surface sub-classes are added, the changes can be made here 
/// without having to touch the Surface class itself.
/// \todo Expand the SurfaceFactory namespace into a singleton class.  Add 
/// pairs of strings and function pointers to an STL map so that the 
/// createSurface methods just do a lookup in the map instead of the clunky
/// if...else construct.  Priority: low.
namespace SurfaceFactory {
  /// Open up the file specified by parameter filename.  The first item in the 
  /// file should be the name of the surface type.  Once that information is 
  /// known, build create and return the appropriate surface.  The client is
  /// responsible to call delete on the Surface* that is returned.
  Surface* createSurface(const std::string& filename);

  /// The parameter type specifies the name of a Surface class.  Instantiate
  /// the appropriate class using parameter surfData, and return the new
  /// surface.  The client is responsible to call delete on the Surface* that
  /// is returned.
  Surface* createSurface(const std::string& type, SurfData* surfData);
}
#endif
