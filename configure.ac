#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

# The version of autoconf must be 2.58 or later
AC_PREREQ([2.58])

# This is Surfpack, version 0.1.  User will be invited to send email
# to the address specified in the third argument when something goes
# wrong.
AC_INIT([Surfpack], [0.1], [mricha@sandia.gov])

# The macros with the AM prefix are instructions to automake.
# Automake will base some of its behavior on things it sees in this file.

# Tells automake that a config.h will be used.  config.h is created by autoconf
# based on config.h.in and the results of the other macros.  As an example,
# suppose one of the macros checks for the BLAS library and tells the C compiler
# to #define BLAS if it exists.  All of those #define results from the macros 
# will appear in a file called config.h.  The file config.h must be #included
# in all of the appropriate source files.  If the AM_CONFIG_HEADER option is not
# used, then all of the #defines will be listed as a command line argument to 
# the compiler for each compilation unit (i.e. g++ -DBLAS).  The results are 
# functionally the same, but it makes it harder to sort through the build
# output when a half page of -D___ arguments appear on the screen for every
# compilation unit.

# The file config.h.in should be generated automatically using the autoheader
# utility.  Autoheader scans configure.ac to see which symbols it may #define.
AM_CONFIG_HEADER([config.h])
AM_INIT_AUTOMAKE

# Check to make sure this file exists as a safety check.  This is merely
# sanity check to make sure you're not running autoconf recklessly
AC_CONFIG_SRCDIR([src/Surface.cpp])

# Add all the necessary m4 directories to aclocal's search path.
AC_CONFIG_MACRO_DIR([m4])

# Checks for one of the following fortran compilers.  Preference is given
# to those listed first.  For example, if the macro is listed as
# 	AC_PROG_F77([xlf f77 g77])
# autoconf will first try to compile a simple fortran program by invoking
# xlf.  If the command succeeds, that means xlf exists and will be used
# as the Fortran 77 compiler.  If the command had failed, autoconf would
# have tried f77.
AC_PROG_F77([xlf f77 g77])

# Checks for a C++ compiler.  Preference is given to the options listed
# first, as in the Fortran example above
echo "CXX=$CXX"
echo "CC=$CC"
AC_PROG_CXX([xlC cxx CC c++ g++])

# This is a kludgy test.  The goal is to turn on the -fpe3 flag for the
# OSF fortran compiler, which allows the program to continue if when  there 
# are floating point errors.  (This behavior is default on many # platforms.)  
# A better idea would be to write a macro that invokes the chosen fortran 
# compiler with this flag and sets FFLAGS to use -fpe3 if the file compiles
# with out any problems.
if test "$F77" = "f77"
then
  if test "$CXX" = "cxx"
  then
    echo "Compiling fortran 77 with f77"
    FFLAGS="-g -O2 -fpe3"
  fi 
fi

# Macro invocation to determine whether references to FORTRAN 77
# subroutines and functions from C/C++ code require underscores.
AC_F77_WRAPPERS

# This macro causes autoconf to tell automake to create Makefiles that
# generate static libraries by default.  The user could still request
# that shared libraries be built by passing an argument to configure:
#  configure --enable-shared.  In other words, this macro only
# sets the DEFAULT.
AC_DISABLE_SHARED

# Checks for the LAPACK library.  This macro was downloaded from
# www.gnu.org/software/ac-archive (current as of August 2004).
# It is included here in acinclude.m4.
ACX_LAPACK

# On some of the compilers I was using, it seemed that even if a directory
# was included in the LD_LIBRARY_PATH, libraries in that directory could 
# not be linked against just by doing a -llib on the link line.  It was
# necessary to also have a -Ldirectory_path argument.  This code expands
# the LD_LIBRARY_PATH into a list of -LDIR1 -LDIR2 linker directives.
# I think I'm just not aware of how it should really be done.  At any rate,
# this code would be more appropriately placed as a macro in acinclude.m4,
# and merely invoked here in configure.ac.
if test -n $LD_LIBRARY_PATH
then 
  WHOLE_LIB_PATH=`echo $LD_LIBRARY_PATH | sed -e 's/^/-L/;s/:/ -L/g'`
else
  WHOLE_LIB_PATH=
fi

# Checks for the flex version and associated library (may not be -lfl)
# that the parser/scanner code needs to link against.  Also check for
# bison.
AC_PROG_LEX
AC_PROG_YACC

# Checks for C++ STL implementation (and ISO standard compliance).
AC_CXX_HAVE_STL
AC_CXX_HAVE_STD

# Macro written by Mark Richards.  Writes a little C++ program and tries
# to link against -lcppunit.
MDR_CXX_CPPUNIT

# The AC_DEFINE macro causes its first argument to always be #defined when
# source code is compiled.
AC_DEFINE([__TESTING_MODE__],[1],
  [Causes additional debug code to be compiled in source.])

# Tells autoconf that libtool will be used
AC_PROG_LIBTOOL

# Sets the flags for the C/C++ preprocessor
CPPFLAGS=""

# Sets the compile line flags for the C++ compiler.  Some of these are
# gcc specific.  Profile-arcs and test-coverage must be activated if you
# want to run the gcov profiler to test code coverage.  The save-temps
# flag tells the compiler to keep the preprocessed source.
#CXXFLAGS="-fprofile-arcs -ftest-coverage -g -save-temps"
CXXFLAGS="$CXXFLAGS -I/usr/local/include -g"

# The cxx compiler on OSF needs the __USE_STD_IOSTREAM flag since
# the sstream classes are used.
if test "$CXX" = "cxx"
then
  echo "Compiling with cxx"
  AC_DEFINE([__USE_STD_IOSTREAM],[],[Define if sstream is used with cxx])
fi

# Autoconf creates a configure script.  When that configure script is run, the
# appropriate Makefile.in files are created.  The .in suffixes are not listed
# explicitly here, but that's what happens.  Autoconf creates the Makefile.in
# files, and then configure creates the Makefiles from the Makefile.in files.
AC_CONFIG_FILES([Makefile src/Makefile])

# Macro written by Mark Richards.  This macro tells autoconf to generate an
# additional Makefile, if and only if the test for CPP Unit succeeded.
MDR_CPPUNIT_MAKEFILE([tests/Makefile])

# For good measure, print out the $LIBS variable at the end.  Any library 
# listed in LIBS is automatically linked against for every executable that is
# created.
echo 'LIBS:' $LIBS

# Mandatory last line.  It tells autoconf to actually generate the appropriate 
# files.
AC_OUTPUT
