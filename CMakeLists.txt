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

# Surfpack requires Boost for serialization
#??? set(Boost_USE_MULTITHREADED OFF)
if(Boost_DIR)
  message("in ${CMAKE_CURRENT_BINARY_DIR} Boost_DIR already set: ${Boost_DIR}")
  include_directories(${Boost_INCLUDE_DIR} ${Boost_INCLUDE_DIRS})
else()
  find_package(Boost 1.37 REQUIRED COMPONENTS serialization)
  if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIR} ${Boost_INCLUDE_DIRS})
  endif()
endif(Boost_DIR)


add_definitions("-D__TESTING_MODE__=1")
remove_definitions("-DHAVE_CONFIG_H")

include_directories(${Surfpack_BINARY_DIR})

FortranCInterface_HEADER(surf77_config.h MACRO_NAMESPACE SURF77_)

install(FILES ${Surfpack_BINARY_DIR}/surf77_config.h DESTINATION include)

set(Surfpack_ENABLE_DEBUG OFF)

#include(CTest)

option(SURFPACK_STANDALONE  "Create a standalone surfpack executable" OFF)
if(SURFPACK_STANDALONE)
  find_package(BISON)
  find_package(FLEX)
endif()

add_subdirectory(packages)
add_subdirectory(src)
#add_subdirectory(tests)
if(SURFPACK_STANDALONE)
  add_subdirectory(interface)
endif()
