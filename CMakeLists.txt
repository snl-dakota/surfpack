cmake_minimum_required(VERSION 2.8)

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

FortranCInterface_HEADER(include/surf77_config.h MACRO_NAMESPACE SURF77_)

set(Surfpack_ENABLE_DEBUG OFF)

#INCLUDE(CTest)

add_subdirectory(packages)
add_subdirectory(src)
#add_subdirectory(tests)

#Keep track of what libraries should be included in libsurfpack
set(SUBLIBS ${SUBLIBS} surfpack)

#Descend into source subdirectories
##SET(SUBLIBS ${SUBLIBS} mars surfaces conmin)
set(SUBLIBS ${SUBLIBS} mars surfaces nkm)

