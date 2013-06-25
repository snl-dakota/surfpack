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

# Needed on MacOSX now that fortran code is split-out into separate library
if(CMAKE_SYSTEM_NAME MATCHES Darwin)
  set(CMAKE_SHARED_LIBRARY_CREATE_Fortran_FLAGS
      "${CMAKE_SHARED_LIBRARY_CREATE_Fortran_FLAGS} -undefined dynamic_lookup")
endif()

option(SURFPACK_STANDALONE  "Create a standalone surfpack executable" OFF)
if(SURFPACK_STANDALONE)
  # By default, require Boost serialization in stand-alone
  if(NOT DEFINED SURFPACK_HAVE_BOOST_SERIALIZATION)
    message(STATUS 
      "SURFPACK_STANDALONE enabling SURFPACK_HAVE_BOOST_SERIALIZATION")
    set(SURFPACK_HAVE_BOOST_SERIALIZATION ON CACHE BOOL 
      "Use Boost serialization library")
  endif()
endif()

option(ENABLE_SPEC_MAINT "Enable Surfpack specification maintenance mode." OFF)
if(ENABLE_SPEC_MAINT AND NOT UNIX)
  message(FATAL_ERROR 
    "Surfpack specification maintenance mode only available on UNIX platforms")
endif()

if(ENABLE_SPEC_MAINT)
  find_package(BISON)
  find_package(FLEX)
endif()

# Surfpack optionally uses Boost for serialization
#??? set(Boost_USE_MULTITHREADED OFF)
option(SURFPACK_HAVE_BOOST_SERIALIZATION "Use Boost serialization library" OFF)
if (SURFPACK_HAVE_BOOST_SERIALIZATION)
  if(Boost_DIR)
    message("in ${CMAKE_CURRENT_BINARY_DIR} Boost_DIR already set: ${Boost_DIR}")
    include_directories(${Boost_INCLUDE_DIR} ${Boost_INCLUDE_DIRS})
  else()
    find_package(Boost 1.37 REQUIRED COMPONENTS serialization)
    if(Boost_FOUND)
      include_directories(${Boost_INCLUDE_DIR} ${Boost_INCLUDE_DIRS})
    endif()
  endif(Boost_DIR)
  add_definitions("-DSURFPACK_HAVE_BOOST_SERIALIZATION")
endif()

add_definitions("-D__TESTING_MODE__=1")
remove_definitions("-DHAVE_CONFIG_H")

include_directories(${Surfpack_BINARY_DIR})

FortranCInterface_HEADER(surf77_config.h MACRO_NAMESPACE SURF77_)

install(FILES ${Surfpack_BINARY_DIR}/surf77_config.h DESTINATION include)

set(Surfpack_ENABLE_DEBUG OFF)

set(Surfpack_NKM_Tests OFF CACHE BOOL "Build the NKM test binaries?")

#include(CTest)

add_subdirectory(packages)
add_subdirectory(src)
#add_subdirectory(tests)
if(SURFPACK_STANDALONE)
  add_subdirectory(interface)
endif()
