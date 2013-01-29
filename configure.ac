#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

# The version of autoconf must be 2.60 or later, due to CONMIN dependency
AC_PREREQ([2.60])

# This is Surfpack, version 1.1.  User will be invited to send email
# to the address specified in the third argument when something goes
# wrong.
AC_INIT([Surfpack],[1.1],[surfpack-developers@development.sandia.gov])

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
AM_INIT_AUTOMAKE([1.9.6 foreign subdir-objects])

# Check to make sure this file exists as a safety check.  This is
# merely a sanity check to make sure you're not running autoconf
# recklessly.
AC_CONFIG_SRCDIR([src/surfpack.cpp])

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

# Surfpack is a C++ code.
AC_LANG([C++])

# Checks for the LAPACK library.  This macro was downloaded from
# http://www.gnu.org/software/ac-archive (current as of August 2004).
# It is included here in the m4 directory.
ACX_BLAS
ACX_LAPACK

# Checks for the flex version and associated library (may not be -lfl)
# that the parser/scanner code needs to link against.  Also check for
# bison.
AC_PROG_LEX
AC_PROG_YACC

# Only build the standalone program if it is specifically requested
AC_ARG_ENABLE(standalone,AS_HELP_STRING([--enable-standalone],[build command-line interpreter (needs lex, yacc)]),[enable_standalone=$enableval],[enable_standalone=no])
AM_CONDITIONAL([BUILD_STANDALONE], [test "x$enable_standalone" = xyes])

# Boost is optional, but don't want to entangle/conditional the probes
AX_BOOST_BASE([1.37], [optional])
if test "x$enable_standalone" = xyes -a "x$succeeded" = "xyes"; then
  AC_DEFINE([SURFPACK_HAVE_BOOST_SERIALIZATION],[1],
            [Whether Surfpack should use Boost serialization library.])
fi

# The AC_DEFINE macro causes its first argument to always be #defined
# when source code is compiled.
AC_DEFINE([__TESTING_MODE__],[1],
  [Causes additional debug code to be compiled in source.])

# Tells autoconf that libtool will be used.
AC_PROG_LIBTOOL

# Checks for library functions.
## GSL not presently needed
#AC_CHECK_LIB([gslcblas],[cblas_ddot])
#AC_CHECK_LIB([gsl],[gsl_ran_fdist_pdf])


# Surfpack-specific package checks. (WJB: NOT needed if we do not use teuchos)
SURFPACK_PACKAGES

# Set up conditional make check target using CppUnit.
AM_PATH_CPPUNIT([1.9.4], [
  AC_MSG_RESULT(

=== Found CppUnit
=== Setting up make check to run tests...)
found_cppunit=yes
],
  AC_MSG_RESULT([

=== Unable to find a recent enough CppUnit
=== Will not run unit testing!!
===
=== To get CppUnit see http://cppunit.sourceforge.net
]))
AM_CONDITIONAL(HAVE_CPPUNIT, test "$found_cppunit" = "yes")

# Autoconf creates a configure script.  When that configure script is
# run, the appropriate Makefile.in files are created.  The .in
# suffixes are not listed explicitly here, but that's what happens.
# Autoconf creates the Makefile.in files, and then configure creates
# the Makefiles from the Makefile.in files.
AC_CONFIG_FILES([Makefile
		 src/Makefile
		 src/interpreter/Makefile
		 src/mars/Makefile
		 src/surfaces/Makefile
		 src/surfaces/nkm/Makefile
		 interface/Makefile
		 packages/Makefile
		 system_tests/Makefile
		 system_tests/instances/Makefile
		 tests/Makefile])

# Mandatory last line.  It tells autoconf to actually generate the
# appropriate files.
AC_OUTPUT
