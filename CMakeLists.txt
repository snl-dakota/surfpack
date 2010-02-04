cmake_minimum_required(VERSION 2.8)

PROJECT("Surfpack" CXX C Fortran)

# WJB: Not sure "Trilinos-style" package/macros are warranted for Surfpack
#INCLUDE(PackageLibraryMacros)
INCLUDE(FortranCInterface)
INCLUDE(CheckIncludeFiles)
# WJB - clueless about this one: INCLUDE(CheckHeaderSTDC)
# WJB - Broader, Compatibility macro might be easier: INCLUDE(TestForSTDNamespace)
INCLUDE(CMakeBackwardCompatibilityCXX)

IF(WIN32)
  FIND_PACKAGE(CLAPACK REQUIRED)
ELSE()
  FIND_PACKAGE(LAPACK REQUIRED)
  #FIND_PACKAGE(BLAS REQUIRED)
ENDIF() 

#INCLUDE(CTest)

#SET( CMAKE_MODULE_PATH
#  ${CMAKE_CURRENT_SOURCE_DIR}/cmake
#  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/utils
#  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/package_arch
##  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/config_tests
#)

#INCLUDE(TrilinosFortranMangling)

#SET(Surfpack_ENABLE_C ON)
SET(Surfpack_ENABLE_CXX ON)
SET(Surfpack_ENABLE_Fortran ON)

CHECK_INCLUDE_FILES(stdlib.h HAVE_STDLIB_H)

INCLUDE_DIRECTORIES("${PROJECT_BINARY_DIR}")

# Package-specific configuration options (WJB: how much "trilinos stuff" is needed?)
#PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)
CONFIGURE_FILE("${PROJECT_SOURCE_DIR}/cmake/surfpack_config.h.in"
  "${PROJECT_BINARY_DIR}/surfpack_config.h"
)

SET(Surfpack_ENABLE_DEBUG OFF)
#SET(Surfpack_ENABLE_EXPLICIT_INSTANTIATION ON)

LIST(APPEND CMAKE_CXX_FLAGS -DHAVE_CONFIG_H)

#SUBDIRS(packages src tests)
SUBDIRS(packages src)

#Keep track of what libraries should be included in libsurfpack
# WJB ToDo: ask gurus about build ordering - would like to build f77 code FIRST!
SET(SUBLIBS surfpack)

#Descend into source subdirectories
SET(SUBLIBS ${SUBLIBS} mars surfaces conmin)

# libsurfpack.la must be built PRIOR to the linkage of test executables
##ADD_LIBRARY(surfpack empty.cxx)
##ADD_LIBRARY(surfpack)
##TARGET_LINK_LIBRARIES(surfpack ${SUBLIBS})

# WJB - DEBUG: look for vars to use for TestForSTDNamespace (cstdlib/std::) and CheckHeaderSTDC (latter is probably not worthy of significant effort)
FILE(WRITE ${CMAKE_CURRENT_BINARY_DIR}/AllVariables.txt "")
GET_CMAKE_PROPERTY(res VARIABLES)
FOREACH(var ${res})
  FILE(APPEND ${CMAKE_CURRENT_BINARY_DIR}/AllVariables.txt
             "${var} \"${${var}}\"\n")
ENDFOREACH(var ${res})

