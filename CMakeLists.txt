cmake_minimum_required(VERSION 2.8)

project("Surfpack" CXX C Fortran)


# Autotools items not currently address in Surfpack's CMake:
# * CPPUnit
# * system_tests


# Set the export name for install targets; parent packages likely want
# to override this to the name of their target
set(ExportTarget ${CMAKE_PROJECT_NAME}-targets CACHE STRING 
  "Name for the export target for ${CMAKE_PROJECT_NAME}")

include(FortranCInterface)

# TODO: Update to use CMake BLAS/LAPACK probes and properly use their output
if(NOT DEFINED BLAS_LIBS OR NOT DEFINED LAPACK_LIBS)
  # Historically on MSVC, tried to use CMake config file approach first.
  # Could probably just use the Unix logic below instead...
  if(MSVC)
    find_package(LAPACK REQUIRED NO_MODULE)
    set(BLAS_LIBS blas)
    set(LAPACK_LIBS lapack)
  else()
    # first check for a system blas and lapack
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
endif() 

# Surfpack optionally depends on Pecos / Teuchos for orthogonal
# matching pursuit solver
if (HAVE_PECOS)
  add_definitions("-DHAVE_PECOS")
  include_directories(${Pecos_BINARY_DIR} ${Pecos_SOURCE_DIR}/util/src
    ${Pecos_SOURCE_DIR}/src)

  # Surfpack Teuchos management; 2 cases to consider, but support same
  # 3 as OPT++ even though this feature can currently only be enabled
  # running under Dakota.
  #  1. Installed Trilinos
  #  2. Teuchos provided elsewhere in build, e.g., Dakota
  #  3. Local Teuchos (separately checked out into trilinos/)

  if (NOT BUILD_IN_TRILINOS)

    # Workaround to skip finding system Trilinos until this probe is
    # simplified and follows find_package semantics. Double negative to
    # preserve historical behavior without overcomplicating things.
    if(NOT DAKOTA_NO_FIND_TRILINOS)
      # First probe for system-installed Trilinos, respecting Trilinos_DIR if set
      find_package(Trilinos QUIET)
    endif()

    if(Trilinos_FOUND)
      if(Trilinos_DIR)
	message(STATUS
	  "Surfpack using previously specified Trilinos in ${Trilinos_DIR}")
      else()
	message(STATUS "Surfpack using external Trilinos")
      endif()
    else()

      # If no parent project configured Teuchos, do so, using Teuchos_DIR if set
      if(Teuchos_DIR)
	message(STATUS
	  "Surfpack using previously specified Teuchos in ${Teuchos_DIR}")
      elseif(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/packages/trilinos)

	set(Teuchos_DIR
          ${CMAKE_CURRENT_BINARY_DIR}/packages/trilinos/packages/teuchos)
	set(Trilinos_ENABLE_Teuchos ON CACHE BOOL
          "Surfpack enabling Trilinos Teuchos" FORCE)

	# Map key Surfpack variables to TriBITS variables
	set( TPL_BLAS_LIBRARIES ${BLAS_LIBS} )
	set( TPL_LAPACK_LIBRARIES ${LAPACK_LIBS} )
	# Surfpack doesn't use any Teuchos MPI features; may want to force off
	#set( TPL_ENABLE_MPI ${DAKOTA_HAVE_MPI} )
	# This doesn't do as name implies; setting OFF doesn't generate Config.cmake
	# at all; doesn't just control whether installed!  Want Config.cmake in build
	#        set(Trilinos_ENABLE_INSTALL_CMAKE_CONFIG_FILES OFF CACHE BOOL
	#	  "Surfpack is the top-level package; don't write Trilinos config files")

	message(STATUS "Surfpack setting Teuchos_DIR to ${Teuchos_DIR}")
	add_subdirectory(packages/trilinos)

      else()
	message(WARNING
	  "Surfpack requires teuchos, but could not find local nor parent copy!")
      endif() # Teuchos_DIR

      # Additional setting to prevent multiple targets with the same name
      set(Trilinos_TARGETS_IMPORTED 1)

      find_package( Teuchos NO_MODULE REQUIRED )

    endif() # Trilinos_DIR

  endif() # NOT BUILD_IN_TRILINOS

  include_directories(${Teuchos_INCLUDE_DIRS})
  link_directories(${Teuchos_LIBRARY_DIRS})

endif()

# WJB - ToDo: introduce teuchos dependency when requested by team
#include_directories(${Surfpack_BINARY_DIR}/packages/teuchos/src
#                    ${Surfpack_SOURCE_DIR}/packages/teuchos/src)


# Needed on MacOSX now that fortran code is split-out into separate library
if(CMAKE_SYSTEM_NAME MATCHES Darwin)
  set(CMAKE_SHARED_LIBRARY_CREATE_Fortran_FLAGS
      "${CMAKE_SHARED_LIBRARY_CREATE_Fortran_FLAGS} -undefined dynamic_lookup")
endif()

# Always need to be on, but Dakota might turn them off
option(HAVE_CONMIN "Build the CONMIN package." ON)
if (NOT HAVE_CONMIN)
  message(SEND_ERROR "Surfpack requires HAVE_CONMIN:BOOL=ON")
endif()
option(HAVE_NCSUOPT "Build the NCSUOPT package." ON)
if (NOT HAVE_NCSUOPT)
  message(SEND_ERROR "Surfpack requires HAVE_NCSUOPT:BOOL=ON")
endif()


option(SURFPACK_ENABLE_SPEC_MAINT 
  "Enable Surfpack specification maintenance mode." OFF)
if(SURFPACK_ENABLE_SPEC_MAINT AND NOT UNIX)
  message(FATAL_ERROR 
    "Surfpack specification maintenance mode only available on UNIX platforms")
endif()

if(SURFPACK_ENABLE_SPEC_MAINT)
  find_package(BISON)
  find_package(FLEX)
endif()

# Surfpack optionally uses Boost for serialization
#??? set(Boost_USE_MULTITHREADED OFF)
option(HAVE_BOOST_SERIALIZATION "Use Boost serialization library" ON)
if (HAVE_BOOST_SERIALIZATION)
  find_package(Boost 1.58 REQUIRED COMPONENTS serialization)
endif()

option(SURFPACK_STANDALONE  "Create a standalone surfpack executable" ON)
if(SURFPACK_STANDALONE AND NOT HAVE_BOOST_SERIALIZATION)
  message(WARNING
    "SURFPACK_STANDALONE is limited without HAVE_BOOST_SERIALIZATION")
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

# Define these libraries at top level so all targets use the same link order

# Core Surfpack libs (linked 1st)
set(SURFPACK_LIBS "surfpack;surfpack_fortran")

# Surfpack packages/ libs (linked 2nd)
set(SURFPACK_TPL_LIBS "ncsuopt;conmin")
# Pecos and Teuchos are most commonly in the build tree
if (HAVE_PECOS)
  list(APPEND SURFPACK_TPL_LIBS ${Pecos_LINK_LIBRARIES} ${Teuchos_LIBRARIES})
endif()

# System libs (linked 3rd)
# Shouldn't need ${FLEX_LIBRARIES} any more due to %noyywrap in lexer.l
set(SURFPACK_SYSTEM_LIBS "${LAPACK_LIBS};${BLAS_LIBS}")


add_subdirectory(src)
#add_subdirectory(tests)
if(SURFPACK_STANDALONE)
  add_subdirectory(interface)
endif()
#add_subdirectory(examples/CInterface)
