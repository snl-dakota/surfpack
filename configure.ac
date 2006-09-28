#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

# The version of autoconf must be 2.58 or later
AC_PREREQ([2.58])

# This is Surfpack, version 1.0.  User will be invited to send email
# to the address specified in the third argument when something goes
# wrong.
AC_INIT([Surfpack],[1.0],[surfpack@scico.sandia.gov])

# The macros with the AM prefix are instructions to automake.
# Automake will base some of its behavior on things it sees in this
# file.

# Tells automake that a surfpack_config.h will be used.  config.h is
# created by autoconf based on surfpack_config.h.in and the results of
# the other macros.  As an example, suppose one of the macros checks
# for the BLAS library and tells the C compiler to #define BLAS if it
# exists.  All of those #define results from the macros will appear in
# a file called surfpack_config.h.  The file config.h must be
# #included in all of the appropriate source files.  If the
# AM_CONFIG_HEADER option is not used, then all of the #defines will
# be listed as a command line argument to the compiler for each
# compilation unit (i.e. g++ -DBLAS).  The results are functionally
# the same, but it makes it harder to sort through the build output
# when a half page of -D___ arguments appear on the screen for every
# compilation unit.

# The file surfpack_config.h.in should be generated automatically
# using the autoheader utility.  Autoheader scans configure.ac to see
# which symbols it may #define.
AM_CONFIG_HEADER([surfpack_config.h])
AM_INIT_AUTOMAKE([1.9 foreign subdir-objects])

# Check to make sure this file exists as a safety check.  This is
# merely a sanity check to make sure you're not running autoconf
# recklessly.
AC_CONFIG_SRCDIR([src/Surface.cpp])

# Add all the necessary m4 directories to aclocal's search path.
AC_CONFIG_MACRO_DIR([m4])

# Checks for one of the following fortran compilers.  Preference is
# given to those listed first.  For example, if the macro is listed as:
# 	AC_PROG_F77([xlf f77 g77])
# autoconf will first try to compile a simple fortran program by
# invoking xlf.  If the command succeeds, that means xlf exists and
# will be used as the Fortran 77 compiler.  If the command had failed,
# autoconf would have tried f77.
#AC_PROG_F77([xlf f77 g77])
AC_PROG_F77

# Checks for a C++ compiler.  Preference is given to the options
# listed first, as in the Fortran example above.
#AC_PROG_CXX([xlC cxx CC c++ g++])
AC_PROG_CXX

# Check required when 'subdir-objects' automake flag is set.
AM_PROG_CC_C_O

# Macro invocation to determine whether references to FORTRAN 77
# subroutines and functions from C/C++ code require underscores.
AC_F77_WRAPPERS

# This macro causes autoconf to tell automake to create Makefiles that
# generate static libraries by default.  The user could still request
# that shared libraries be built by passing an argument to configure:
#      configure --enable-shared
# In other words, this macro only sets the DEFAULT.
AC_DISABLE_SHARED

# Checks for the LAPACK library.  This macro was downloaded from
# http://www.gnu.org/software/ac-archive (current as of August 2004).
# It is included here in the m4 directory.
ACX_LAPACK
AM_CONDITIONAL([NEEDS_EXTRA_LAPACK_ROUTINES],[test "x$acx_lapack_ok" != xyes])

# Checks for the flex version and associated library (may not be -lfl)
# that the parser/scanner code needs to link against.  Also check for
# bison.
AC_PROG_LEX
AC_PROG_YACC

# If some kind of lex and yacc are found, build the standalone executable 
ACX_BUILD_STANDALONE

# Checks for C++ ISO standard compliance.
AC_CXX_HAVE_STD

# The AC_DEFINE macro causes its first argument to always be #defined
# when source code is compiled.
AC_DEFINE([__TESTING_MODE__],[1],
  [Causes additional debug code to be compiled in source.])

# Tells autoconf that libtool will be used.
AC_PROG_LIBTOOL

# The cxx compiler on OSF needs the __USE_STD_IOSTREAM flag since the
# sstream classes are used.  We also need local template instantiation.
AC_CXX_CPPFLAGS_STD_LANG([CXXFLAGS])
case "$ac_cv_cxx_compiler_vendor" in
    compaq)
	AX_CXXFLAGS_OSF_OPTION([-tlocal])
    ;;
    *)
    ;;
esac

# Set up conditional make check target using CppUnit.
AM_PATH_CPPUNIT([1.9.4], [
  AC_MSG_RESULT(

*** Found CppUnit
*** Setting up make check to run tests...)
found_cppunit=yes
],
  AC_MSG_RESULT([

*** Unable to find a recent enough CppUnit
*** Will not run unit testing!!
***
*** To get CppUnit see http://cppunit.sourceforge.net
]))
AM_CONDITIONAL(HAVE_CPPUNIT, test "$found_cppunit" = "yes")


# Checks for library functions.
AC_CHECK_FUNCS([gettimeofday drand48 srand48])

# Autoconf creates a configure script.  When that configure script is
# run, the appropriate Makefile.in files are created.  The .in
# suffixes are not listed explicitly here, but that's what happens.
# Autoconf creates the Makefile.in files, and then configure creates
# the Makefiles from the Makefile.in files.
AC_CONFIG_FILES([Makefile
		 src/Makefile
		 src/ann/Makefile
		 src/conmin/Makefile
		 src/interpreter/Makefile
		 src/kriging/Makefile
		 src/lapack/Makefile
		 src/mars/Makefile
		 src/surfaces/Makefile
		 interface/Makefile
		 system_tests/Makefile
		 system_tests/instances/Makefile
		 tests/Makefile])

# Mandatory last line.  It tells autoconf to actually generate the
# appropriate files.
AC_OUTPUT
