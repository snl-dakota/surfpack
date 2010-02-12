cmake_minimum_required(VERSION 2.8)

PROJECT("Surfpack" CXX C Fortran)

# WJB - per BMA review, "Trilinos-style pkging is history: INCLUDE(PackageLibraryMacros)
INCLUDE(FortranCInterface)
INCLUDE(CheckIncludeFiles)
####INCLUDE(CMakeBackwardCompatibilityCXX)

IF(WIN32)
  FIND_PACKAGE(CLAPACK REQUIRED)
ELSE()
  FIND_PACKAGE(LAPACK REQUIRED)
ENDIF() 

#INCLUDE(CTest)

SET(Surfpack_ENABLE_CXX ON)
SET(Surfpack_ENABLE_Fortran ON)

INCLUDE_DIRECTORIES("${PROJECT_BINARY_DIR}")

# WJB: hope to eradicate autotools-specific build techiques... (maybe someday)
CONFIGURE_FILE("${PROJECT_SOURCE_DIR}/cmake/surfpack_config.h.in"
  "${PROJECT_BINARY_DIR}/surfpack_config.h"
)

SET(Surfpack_ENABLE_DEBUG OFF)
# WJB: seems a bit Trilinos-specific, but will defer...
#SET(Surfpack_ENABLE_EXPLICIT_INSTANTIATION ON)

# WJB: yaHope to eradicate autotools-specific build techiques... (maybe someday)
LIST(APPEND CMAKE_CXX_FLAGS -DHAVE_CONFIG_H)

#SUBDIRS(packages src tests)
SUBDIRS(packages src)

#Keep track of what libraries should be included in libsurfpack
# WJB ToDo: ask gurus about build ordering - would like to build f77 code FIRST!
SET(SUBLIBS surfpack)

#Descend into source subdirectories
SET(SUBLIBS ${SUBLIBS} mars surfaces conmin)

