
include_directories(${Surfpack_BINARY_DIR}
                    ${Surfpack_SOURCE_DIR}/src
                    ${Surfpack_SOURCE_DIR}/src/surfaces)
# WJB - ToDo: introduce interpreter dependency when I am ready for LEX/YACC
#                    ${Surfpack_SOURCE_DIR}/src/interpreter
#
# WJB - ToDo: introduce teuchos dependency when requested by team
#include_directories(${Surfpack_BINARY_DIR}/packages/teuchos/src
#                    ${Surfpack_SOURCE_DIR}/packages/teuchos/src)

add_subdirectory(mars)
add_subdirectory(surfaces)

## Surfpack library sources.
set(Surfpack_src
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

# WJB - ToDo: look into CMake template support features: SurfpackMatrix.tpl

add_library(surfpack ${Surfpack_src})
#target_link_libraries(surfpack_src mars teuchos conmin ${LAPACK_LIBRARIES})
install(TARGETS surfpack EXPORT ${ExportTarget} DESTINATION lib)

install(FILES SurfpackMatrix.h surfpack_system_headers.h DESTINATION include)
