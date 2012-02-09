cmake_minimum_required(VERSION 2.8)

project("Surfpack" CXX C Fortran)

# Set the export name for install targets; parent packages likely want
# to override this to the name of their target
set(ExportTarget ${CMAKE_PROJECT_NAME}-targets CACHE STRING 
  "Name for the export target for ${CMAKE_PROJECT_NAME}")

include(FortranCInterface)

# first check for a system blas and lapack
if(NOT DEFINED BLAS_LIBS OR NOT DEFINED LAPACK_LIBS)

  if(NOT DEFINED BLAS_LIBS)
    find_library(BLAS_LIBS blas)
  endif()
  if(NOT DEFINED LAPACK_LIBS)
    find_library(LAPACK_LIBS lapack)
  endif()
  if(NOT BLAS_LIBS OR NOT LAPACK_LIBS)
    # if not a system blas and lapack, then look for a cmake built LAPACK
    # with find_package
    find_package(LAPACK REQUIRED NO_MODULE)
    set(BLAS_LIBS blas)
    set(LAPACK_LIBS lapack)
  endif()

endif() 

add_definitions("-D__TESTING_MODE__=1")
remove_definitions("-DHAVE_CONFIG_H")

include_directories(${Surfpack_BINARY_DIR})

FortranCInterface_HEADER(surf77_config.h MACRO_NAMESPACE SURF77_)

install(FILES ${Surfpack_BINARY_DIR}/surf77_config.h DESTINATION include)

set(Surfpack_ENABLE_DEBUG OFF)

#include(CTest)

add_subdirectory(packages)
add_subdirectory(src)
#add_subdirectory(tests)

