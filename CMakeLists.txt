cmake_minimum_required(VERSION 2.8)

PROJECT("Surfpack" CXX C Fortran)

INCLUDE(FortranCInterface)
INCLUDE(CheckIncludeFiles)

IF(WIN32)
  FIND_PACKAGE(CLAPACK REQUIRED)
ELSE()
  FIND_PACKAGE(LAPACK REQUIRED)
ENDIF() 

SET(Surfpack_ENABLE_CXX ON)
SET(Surfpack_ENABLE_Fortran ON)

INCLUDE_DIRECTORIES("${PROJECT_BINARY_DIR}")

CONFIGURE_FILE("${PROJECT_SOURCE_DIR}/cmake/surfpack_config.h.in"
  "${PROJECT_BINARY_DIR}/surfpack_config.h"
)

SET(Surfpack_ENABLE_DEBUG OFF)

LIST(APPEND CMAKE_CXX_FLAGS -DHAVE_CONFIG_H)

#INCLUDE(CTest)

#SUBDIRS(packages src tests)
SUBDIRS(packages src)

#Keep track of what libraries should be included in libsurfpack
# WJB ToDo: ask gurus about build ordering - would like to build f77 code FIRST!
SET(SUBLIBS surfpack)

#Descend into source subdirectories
SET(SUBLIBS ${SUBLIBS} mars surfaces conmin)

