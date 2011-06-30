
include_directories(${Surfpack_BINARY_DIR}
                    ${Surfpack_SOURCE_DIR}/src)

# WJB - ToDo: introduce interpreter dependency when I am ready for LEX/YACC
#                    ${Surfpack_SOURCE_DIR}/src/interpreter
#
# WJB - ToDo: introduce teuchos dependency when requested by team
#include_directories(${Surfpack_BINARY_DIR}/packages/teuchos/src
#                    ${Surfpack_SOURCE_DIR}/packages/teuchos/src)

set(local_library surfpack)

## Top-level Surfpack library sources.
set(lib${local_library}_sources
   AxesBounds.cpp
   AxesBounds.h
   MersenneTwister.h
   ModelFactory.cpp
   ModelFactory.h
   SurfData.cpp
   SurfData.h
   surfpack.cpp
   surfpack.h
   SurfpackMatrix.h
   surfpack_system_headers.h
   SurfpackInterface.h
   SurfpackInterface.cpp
   SurfpackParserArgs.h
   SurfpackParserArgs.cpp
   SurfPoint.cpp
   SurfPoint.h
   Conmin.cpp
   Conmin.h
   ModelFitness.cpp
   ModelFitness.h
)

# Build a single libsurfpack from all subdirectories
set(local_src_subdirs mars surfaces surfaces/nkm)

# Aggregate all the sources into a single library for installation
foreach(sd_it ${local_src_subdirs})

  add_subdirectory(${sd_it})

  # grab the trailing part of any paths
  string(REGEX MATCH "[^/]+$" last_path_name ${sd_it})

  # extract the variables ${last_path_name}_sources from the subdir's env
  get_directory_property(sd_sources DIRECTORY ${sd_it} 
    DEFINITION ${last_path_name}_sources)

  # add the sources to the library's list and includes
  foreach ( src_file ${sd_sources} )
    list ( APPEND lib${local_library}_sources "${sd_it}/${src_file}" )
  endforeach()

  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/${sd_it})

endforeach()

# WJB - ToDo: look into CMake template support features: SurfpackMatrix.tpl

add_library(${local_library} ${lib${local_library}_sources})
#target_link_libraries(surfpack_src mars teuchos conmin ${LAPACK_LIBRARIES})
install(TARGETS ${local_library} EXPORT ${ExportTarget} DESTINATION lib)

install(FILES SurfpackMatrix.h surfpack_system_headers.h DESTINATION include)
