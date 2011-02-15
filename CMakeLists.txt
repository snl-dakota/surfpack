cmake_minimum_required(VERSION 2.6)

project("Surfpack" CXX C Fortran)

include(FortranCInterface)

if(WIN32)
  find_package(CLAPACK REQUIRED)
else()
  # first check for a system blas and lapack
  find_library(BLAS_LIBS blas)
  find_library(LAPACK_LIBS lapack)
  if(NOT BLAS_LIBS OR NOT LAPACK_LIBS)
    # if not a system blas and lapack, then
    # look for a cmake built LAPACK with
    # find_package
    find_package(LAPACK REQUIRED NO_MODULE)
    set(BLAS_LIBS blas)
    set(LAPACK_LIBS lapack)
  endif()
#  FIND_PACKAGE(BLAS REQUIRED)
#  FIND_PACKAGE(LAPACK REQUIRED)
endif() 

include_directories("${PROJECT_BINARY_DIR}")

configure_file("${PROJECT_SOURCE_DIR}/cmake/surfpack_config.h.in"
  "${PROJECT_BINARY_DIR}/surfpack_config.h"
  @ONLY
)

add_definitions("-DHAVE_CONFIG_H")
#add_definitions("-DF77_FUNC_=@FortranCInterface_GLOBAL_MACRO@")

#SET(Surfpack_ENABLE_CXX ON)
#SET(Surfpack_ENABLE_Fortran ON)
set(Surfpack_ENABLE_DEBUG OFF)

#INCLUDE(CTest)

if ( NOT CONMIN_DIR )
  set( CONMIN_DIR ${CMAKE_CURRENT_BINARY_DIR}/packages/CONMIN )
  message( "Setting CONMIN_DIR to ${CONMIN_DIR}" )
  set(SUBLIBS ${SUBLIBS} conmin)
else()
  message( "in ${CMAKE_CURRENT_BINARY_DIR} CONMIN_DIR already set to ${CONMIN_DIR}" )
endif()

add_subdirectory(packages)
add_subdirectory(src)
#add_subdirectory(tests)

#Keep track of what libraries should be included in libsurfpack
# WJB ToDo: ask gurus about build ordering - would like to build f77 code FIRST!
set(SUBLIBS ${SUBLIBS} surfpack)

#Descend into source subdirectories
##SET(SUBLIBS ${SUBLIBS} mars surfaces conmin)
set(SUBLIBS ${SUBLIBS} mars surfaces nkm)
