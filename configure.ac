#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.57)
AC_INIT(Surfpack++, 0.1, mricha@sandia.gov)
AM_INIT_AUTOMAKE
# Check to make sure this file exists as a safety check
AC_CONFIG_SRCDIR([working/snapshot.m])

# Checks for programs.
AC_PROG_AWK
AC_PROG_CXX
AC_PROG_F77
AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_YACC
AC_PROG_MAKE_SET
AC_PROG_LIBTOOL

# Checks for libraries.
## FIXME: Replace `main' with a function in `-lF':
#AC_CHECK_LIB([F], [main])
## FIXME: Replace `main' with a function in `-lICE':
#AC_CHECK_LIB([ICE], [main])
## FIXME: Replace `main' with a function in `-lSM':
#AC_CHECK_LIB([SM], [main])
## FIXME: Replace `main' with a function in `-lX11':
#AC_CHECK_LIB([X11], [main])
## FIXME: Replace `main' with a function in `-lXext':
#AC_CHECK_LIB([Xext], [main])
## FIXME: Replace `main' with a function in `-lXm':
#AC_CHECK_LIB([Xm], [main])
## FIXME: Replace `main' with a function in `-lXmu':
#AC_CHECK_LIB([Xmu], [main])
## FIXME: Replace `main' with a function in `-lXp':
#AC_CHECK_LIB([Xp], [main])
## FIXME: Replace `main' with a function in `-lXpm':
#AC_CHECK_LIB([Xpm], [main])
## FIXME: Replace `main' with a function in `-lXt':
#AC_CHECK_LIB([Xt], [main])
## FIXME: Replace `main' with a function in `-lann':
#AC_CHECK_LIB([ann], [main])
## FIXME: Replace `main' with a function in `-lblas':
#AC_CHECK_LIB([blas], [main])
## FIXME: Replace `main' with a function in `-lconmin':
#AC_CHECK_LIB([conmin], [main])
## FIXME: Replace `main' with a function in `-ldakota':
#AC_CHECK_LIB([dakota], [main])
## FIXME: Replace `main' with a function in `-ldl':
#AC_CHECK_LIB([dl], [main])
## FIXME: Replace `main' with a function in `-lfl':
#AC_CHECK_LIB([fl], [main])
## FIXME: Replace `main' with a function in `-lg2c':
#AC_CHECK_LIB([g2c], [main])
## FIXME: Replace `main' with a function in `-lkriging':
#AC_CHECK_LIB([kriging], [main])
## FIXME: Replace `main' with a function in `-llapack':
#AC_CHECK_LIB([lapack], [main])
## FIXME: Replace `main' with a function in `-lm':
#AC_CHECK_LIB([m], [main])
## FIXME: Replace `main' with a function in `-lmars':
#AC_CHECK_LIB([mars], [main])
## FIXME: Replace `main' with a function in `-lpgc':
#AC_CHECK_LIB([pgc], [main])
## FIXME: Replace `main' with a function in `-lpgf90':
#AC_CHECK_LIB([pgf90], [main])
## FIXME: Replace `main' with a function in `-lpgf902':
#AC_CHECK_LIB([pgf902], [main])
## FIXME: Replace `main' with a function in `-lpgf90_rpm1':
#AC_CHECK_LIB([pgf90_rpm1], [main])
## FIXME: Replace `main' with a function in `-lpgftnrtl':
#AC_CHECK_LIB([pgftnrtl], [main])
## FIXME: Replace `main' with a function in `-lpmpich':
#AC_CHECK_LIB([pmpich], [main])
## FIXME: Replace `main' with a function in `-lpthread':
#AC_CHECK_LIB([pthread], [main])
## FIXME: Replace `main' with a function in `-lsurfpack':
#AC_CHECK_LIB([surfpack], [main])

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST

# Checks for library functions.
AC_FUNC_MALLOC
AC_HEADER_STDC
AC_CHECK_FUNCS([floor pow sqrt])

AC_CONFIG_FILES([Makefile
		 src/Makefile
		 src/conmin/Makefile
		 src/kriging/Makefile
		 src/mars/Makefile
                 src/ann/Makefile
		 interface/Makefile
		 tests/Makefile])
AC_DEFINE([__TESTING_MODE__],[])
AC_OUTPUT
