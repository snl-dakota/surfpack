#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.58)
AC_INIT(Surfpack, 0.1, mricha@sandia.gov)
AM_CONFIG_HEADER(config.h)
AM_INIT_AUTOMAKE
# Check to make sure this file exists as a safety check
AC_CONFIG_SRCDIR([src/Surface.cpp])

# Checks for programs.
#AC_PROG_AWK
AC_PROG_CXX([CC cxx xlC g++ c++ gcc])
#AC_CXX_NAMESPACES
#AC_CXX_HAVE_STL
#AC_CXX_BOOL
#AC_CXX_MUTABLE

#AC_DISABLE_SHARED

AC_LANG_CPLUSPLUS
ACX_LAPACK
WHOLE_LIB_PATH=`echo $LD_LIBRARY_PATH | sed -e 's/^/-L/;s/:/ -L/g'`
#Check for linear algebra libraries needed on the AIX
#AC_CHECK_LIB([essl], [dgels_])
#
##Check for the linear algebra libraries that are needed on IRIX
#AC_CHECK_LIB([scs],[dgels_])
#
##Chech for the linear algebra libraries that should be on OSF
#AC_CHECK_LIB([dxml],[dgels_])
AC_CHECK_LIB([cxxstd],[strcpy])
#
##Check for the linear algebra libraries that are needed on Sun/Solaris
AC_CHECK_LIB([sunmath],[logf],[],[],[$WHOLE_LIB_PATH])
AC_CHECK_LIB([F77],[etime_],[],[],[$WHOLE_LIB_PATH])
AC_CHECK_LIB([fsu],[__c_exp],[],[],[$WHOLE_LIB_PATH])
AC_CHECK_LIB([sunperf],[dgels_],[],[],[$WHOLE_LIB_PATH])

#
##Check for the linear algebra libraries on Linux machines
#AC_CHECK_LIB([m],[atanh],[],[],[$WHOLE_LIB_PATH])
#AC_CHECK_LIB([g2c],[c_cos],[],[],[$WHOLE_LIB_PATH])
#AC_CHECK_LIB([blas],[ddot_],[],[],[$WHOLE_LIB_PATH])
#AC_CHECK_LIB([lapack],[dgels_],[],[],[$WHOLE_LIB_PATH])
##AC_CHECK_LIB([gslcblas],[main],[],[],[$WHOLE_LIB_PATH])
##AC_CHECK_LIB([gsl],[main],[],[],[$WHOLE_LIB_PATH])

#Check for the linear algebra libraries that are needed on IRIX
AC_CHECK_LIB([scs],[dgels_])

#Chech for the linear algebra libraries that should be on OSF
AC_CHECK_LIB([dxml],[dgels_])


#Check for libraries needed for CPP Unit test suite
AC_CHECK_LIB([dl],[dlclose],[],[],[$WHOLE_LIB_PATH])
include_tests=yes
AC_CHECK_LIB([cppunit],[main],[],[include_tests=no],[$WHOLE_LIB_PATH])
AM_CONDITIONAL(INCLUDE_TESTS, test $include_tests = yes)

AC_DEFINE([__TESTING_MODE__],[],[1])


AC_PROG_F77
#AC_PROG_CC
#AC_PROG_INSTALL
#AC_PROG_LN_S
#AC_PROG_YACC
#AC_PROG_MAKE_SET
AC_PROG_LIBTOOL

CPPFLAGS="-DDATADIR='\"\$(top_srcdir)/share\"' "
#CXXFLAGS="-fprofile-arcs -ftest-coverage -g -save-temps"
CXXFLAGS="-g "
if test $CXX = "cxx"
then
  echo "Compiling with cxx"
  CXXFLAGS="-g -std strict_ansi"
fi

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
#AC_HEADER_STDBOOL
#AC_C_CONST

# Checks for library functions.
#AC_FUNC_MALLOC
#AC_HEADER_STDC
#AC_CHECK_FUNCS([floor pow sqrt])
AC_CONFIG_FILES([Makefile
		 src/Makefile
		 src/conmin/Makefile
		 src/kriging/Makefile
		 src/mars/Makefile
                 src/ann/Makefile
		 interface/Makefile])
if test $include_tests = yes; then
  AC_CONFIG_FILES([tests/Makefile])
fi
echo 'LIBS:' $LIBS
AC_OUTPUT
