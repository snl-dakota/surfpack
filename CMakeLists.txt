cmake_minimum_required(VERSION 2.6)

PROJECT("Surfpack" CXX C Fortran)

INCLUDE(FortranCInterface)

IF(WIN32)
  FIND_PACKAGE(CLAPACK REQUIRED)
ELSE()
  FIND_PACKAGE(BLAS REQUIRED)
  FIND_PACKAGE(LAPACK REQUIRED)
ENDIF() 

INCLUDE_DIRECTORIES("${PROJECT_BINARY_DIR}")

CONFIGURE_FILE("${PROJECT_SOURCE_DIR}/cmake/surfpack_config.h.in"
  "${PROJECT_BINARY_DIR}/surfpack_config.h"
  @ONLY
)

ADD_DEFINITIONS("-DHAVE_CONFIG_H")

#SET(Surfpack_ENABLE_CXX ON)
#SET(Surfpack_ENABLE_Fortran ON)
SET(Surfpack_ENABLE_DEBUG OFF)

#INCLUDE(CTest)

#SUBDIRS(packages src tests)
SUBDIRS(packages src)

#Keep track of what libraries should be included in libsurfpack
# WJB ToDo: ask gurus about build ordering - would like to build f77 code FIRST!
SET(SUBLIBS surfpack)

#Descend into source subdirectories
SET(SUBLIBS ${SUBLIBS} mars surfaces conmin)

