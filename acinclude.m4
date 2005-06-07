dnl Available from the GNU Autoconf Macro Archive at:
dnl http://www.gnu.org/software/ac-archive/htmldoc/acx_blas.html
dnl

dnl The macros in this file (acinclude.m4) are included in the aclocal.m4 file
dnl that is generated when aclocal is run.  The macros can then be invoked from
dnl the configure.ac file in this directory

dnl This macro checks for the presence of CPP Unit.  Specifically, it writes a
dnl short program that #includes a CPP Unit header file and declares a class
dnl that derives from a CPP Unit base class.  If that program is successfully 
dnl compiled and then successfully linked against the libraries -lcppunit and 
dnl -ldl, then a HAVE_CPPUNIT variable will be #defined, as well as an automake
dnl conditional INCLUDE_TESTS.  The variable CPPUNIT_LIBS will contain the list
dnl of libraries that need to be linked in with any code that uses CPP Unit
dnl classes or functions (i.e., -lcppunit and -ldl).  The variable can be
dnl dereferenced in Makefile.am files as $(CPPUNIT_LIBS).
AC_DEFUN([MDR_CXX_CPPUNIT], [
  AC_PREREQ(2.58)
  AC_CACHE_CHECK([for CPPUnit libraries], ac_cv_cxx_cppunit,
    AC_LANG_SAVE
    mdr_save_CXXFLAGS=$CXXFLAGS
    CXXFLAGS="-L/usr/netpub/cppunit/lib -I/usr/netpub/cppunit/include -lcppunit -ldl"
    dnl CXXFLAGS="-I/usr/netpub/cppunit/include -lcppunit -ldl"
    AC_LANG(C++)
    AC_LINK_IFELSE(
      [AC_LANG_PROGRAM(
      	[[#include <cppunit/extensions/HelperMacros.h>
      
      	 class MyTest : public CppUnit::TestCase
      	 {
      	   MyTest(std::string name) : CppUnit::TestCase( name ) {}
      	   void runTest() { CPPUNIT_ASSERT( 1 == 1); }
      	 };]],
      	[[//no body]])],
      [ac_cv_cxx_cppunit=yes], [ac_cv_cxx_cppunit=no])
    AC_LANG_RESTORE
    CXXFLAGS="$mdr_save_CXXFLAGS")
  if test "$ac_cv_cxx_cppunit" = yes
  then
    CPPUNIT_LIBS="-lcppunit -ldl"
    AC_SUBST([CPPUNIT_LIBS])
    AC_DEFINE(HAVE_CPPUNIT,,[
      define if CPPUnit libraries are available for linking])
  fi
  AM_CONDITIONAL(INCLUDE_TESTS, test "$ac_cv_cxx_cppunit" = yes)])	

dnl This macro will add its first argument ($1) to the list of files that the 
dnl configure script should generate, if and only if the test for CPP Unit
dnl libraries passed.
AC_DEFUN([MDR_CPPUNIT_MAKEFILE], [
  AC_CACHE_CHECK(whether CPPUnit tests should be included, 
    ac_cv_cxx_cppunit, [])
  if test "$ac_cv_cxx_cppunit" = yes
  then
       AC_CONFIG_FILES([$1])
  fi])

dnl If C++ programs invoke fortran routines, the fortran routines must be 
dnl declared extern "C" in the C++ code (to disable name-mangling). Depending
dnl on the platform an underscore may need to be appended to the name of the
dnl fortran calls.  This macro will cause C2F77_CALLS_NEED_UNDERSCORE to be
dnl #defined iff the trailing underscore is required.
AC_DEFUN([MDR_C2F77_UNDERSCORE], [
  AC_CACHE_CHECK(
    whether calls to fortran from C++ require trailing underscore, 
    ac_cv_add_underscore, [
      AC_F77_FUNC([func1])
      if test "$func1" = "func1_" 
      then
        ac_cv_add_underscore=yes
      else
        ac_cv_add_underscore=no
      fi]
  )
  if test "$ac_cv_add_underscore" = "yes"
  then 
      AC_DEFINE([C2F77_CALLS_NEED_UNDERSCORE],[],[
        "Define if function calls from C++ to F77 need a trailing underscore"])
  fi
  ])

dnl Checks for a library with BLAS routines.  This macro was downloaded from
dnl www.gnu.org/software/ac-archive (current August 2004).  Some modifications
dnl were made for the Solaris platform.  This macro is inherently kludgy, 
dnl because it requires a knowledge of what the names of the BLAS libraries
dnl are on each platform of interest.  The macro will not work for any new
dnl platforms, unless the BLAS librar(y/ies) on that platform happen to have
dnl the same names as the libraries on some other platform.  Still, this is
dnl probably the best that autoconf can be expected to do.
AC_DEFUN([ACX_BLAS], [
  AC_PREREQ(2.50)
  AC_REQUIRE([AC_F77_LIBRARY_LDFLAGS])
  acx_blas_ok=no
  
  AC_ARG_WITH(blas,
  	[AC_HELP_STRING([--with-blas=<lib>], [use BLAS library <lib>])])
  case $with_blas in
  	yes | "") ;;
  	no) acx_blas_ok=disable ;;
  	-* | */* | *.a | *.so | *.so.* | *.o) BLAS_LIBS="$with_blas" ;;
  	*) BLAS_LIBS="-l$with_blas" ;;
  esac
  
  # Get fortran linker names of BLAS functions to check for.
  AC_F77_FUNC(sgemm)
  AC_F77_FUNC(dgemm)
  
  acx_blas_save_LIBS="$LIBS"
  LIBS="$LIBS $FLIBS"
  
  # First, check BLAS_LIBS environment variable
  if test $acx_blas_ok = no; then
  if test "x$BLAS_LIBS" != x; then
  	save_LIBS="$LIBS"; LIBS="$BLAS_LIBS $LIBS"
  	AC_MSG_CHECKING([for $sgemm in $BLAS_LIBS])
  	AC_TRY_LINK_FUNC($sgemm, [acx_blas_ok=yes], [BLAS_LIBS=""])
  	AC_MSG_RESULT($acx_blas_ok)
  	LIBS="$save_LIBS"
  fi
  fi
  
  # BLAS linked to by default?  (happens on some supercomputers)
  if test $acx_blas_ok = no; then
  	save_LIBS="$LIBS"; LIBS="$LIBS"
  	AC_CHECK_FUNC($sgemm, [acx_blas_ok=yes])
  	LIBS="$save_LIBS"
  fi
  
  # BLAS in ATLAS library? (http://math-atlas.sourceforge.net/)
  if test $acx_blas_ok = no; then
  	AC_CHECK_LIB(atlas, ATL_xerbla,
  		[AC_CHECK_LIB(f77blas, $sgemm,
  		[AC_CHECK_LIB(cblas, cblas_dgemm,
  			[acx_blas_ok=yes
  			 BLAS_LIBS="-lcblas -lf77blas -latlas"],
  			[], [-lf77blas -latlas])],
  			[], [-latlas])])
  fi
  
  # BLAS in PhiPACK libraries? (requires generic BLAS lib, too)
  if test $acx_blas_ok = no; then
  	AC_CHECK_LIB(blas, $sgemm,
  		[AC_CHECK_LIB(dgemm, $dgemm,
  		[AC_CHECK_LIB(sgemm, $sgemm,
  			[acx_blas_ok=yes; BLAS_LIBS="-lsgemm -ldgemm -lblas"],
  			[], [-lblas])],
  			[], [-lblas])])
  fi
  
  # BLAS in Alpha CXML library?
  if test $acx_blas_ok = no; then
  	AC_CHECK_LIB(cxml, $sgemm, [acx_blas_ok=yes;BLAS_LIBS="-lcxml"])
  fi
  
  # BLAS in Alpha DXML library? (now called CXML, see above)
  if test $acx_blas_ok = no; then
  	AC_CHECK_LIB(dxml, $sgemm, [acx_blas_ok=yes;BLAS_LIBS="-ldxml"])
  fi
  
  # BLAS in Sun Performance library?
  if test $acx_blas_ok = no
  then
      AC_CHECK_LIB([sunmath],[logf],[
        AC_CHECK_LIB([F77],[etime_],[
          AC_CHECK_LIB([fsu],[__c_exp],[
            AC_CHECK_LIB([sunperf],[dgels_],[
              LAPACK_LIBS="-lsunperf -lfsu -lF77 -lsunmath";acx_blas_ok=yes
  	  ],[],[-lfsu -lF77 -lsunmath $WHOLE_LIB_PATH])
  	],[],[-lF77 -lsunmath $WHOLE_LIB_PATH])
        ],[],[-lsunmath $WHOLE_LIB_PATH])
      ],[],[$WHOLE_LIB_PATH])
  fi
  
  # BLAS in SCSL library?  (SGI/Cray Scientific Library)
  if test $acx_blas_ok = no; then
  	AC_CHECK_LIB(scs, $sgemm, [acx_blas_ok=yes; BLAS_LIBS="-lscs"])
  fi
  
  # BLAS in SGIMATH library?
  if test $acx_blas_ok = no; then
  	AC_CHECK_LIB(complib.sgimath, $sgemm,
  		     [acx_blas_ok=yes; BLAS_LIBS="-lcomplib.sgimath"])
  fi
  
  # BLAS in IBM ESSL library? (requires generic BLAS lib, too)
  if test $acx_blas_ok = no; then
  	AC_CHECK_LIB(blas, $sgemm,
  		[AC_CHECK_LIB(essl, $sgemm,
  			[mdr_essl=yes; echo "mdr_essl = yes";
  			 acx_blas_ok=yes; BLAS_LIBS="-lessl -lblas"],
  			[mdr_essl=no; echo "mdr_essl = no"], [-lblas $FLIBS])])
  fi
  
  # Generic BLAS library?
  if test $acx_blas_ok = no; then
  	AC_CHECK_LIB(blas, $sgemm, [acx_blas_ok=yes; BLAS_LIBS="-lblas"])
  fi
  
  AM_CONDITIONAL([NEEDS_EXTRA_LAPACK_ROUTINES], [test "$mdr_essl" = "yes"])
  AC_SUBST(BLAS_LIBS)
  
  LIBS="$acx_blas_save_LIBS"
  
  # Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
  if test x"$acx_blas_ok" = xyes; then
          ifelse([$1],,AC_DEFINE(HAVE_BLAS,1,[Define if you have a BLAS library.]),[$1])
          :
  else
          acx_blas_ok=no
          $2
  fi
  ])dnl ACX_BLAS
dnl Available from the GNU Autoconf Macro Archive at:
dnl http://www.gnu.org/software/ac-archive/htmldoc/acx_lapack.html
dnl

dnl Checks for LAPACK library routines.  This is also an inherently kludgy 
dnl macro.  See comments about the BLAS macro above.
AC_DEFUN([ACX_LAPACK], [
  AC_REQUIRE([ACX_BLAS])
  acx_lapack_ok=no
  
  AC_ARG_WITH(lapack,
          [AC_HELP_STRING([--with-lapack=<lib>], [use LAPACK library <lib>])])
  case $with_lapack in
          yes | "") ;;
          no) acx_lapack_ok=disable ;;
          -* | */* | *.a | *.so | *.so.* | *.o) LAPACK_LIBS="$with_lapack" ;;
          *) LAPACK_LIBS="-l$with_lapack" ;;
  esac
  
  # Get fortran linker name of LAPACK function to check for.
  AC_F77_FUNC(cheev)
  
  # We cannot use LAPACK if BLAS is not found
  if test "x$acx_blas_ok" != xyes; then
          acx_lapack_ok=noblas
  fi
  
  # First, check LAPACK_LIBS environment variable
  if test "x$LAPACK_LIBS" != x; then
          save_LIBS="$LIBS"; LIBS="$LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS"
          AC_MSG_CHECKING([for $cheev in $LAPACK_LIBS])
          AC_TRY_LINK_FUNC($cheev, [acx_lapack_ok=yes], [LAPACK_LIBS=""])
          AC_MSG_RESULT($acx_lapack_ok)
          LIBS="$save_LIBS"
          if test acx_lapack_ok = no; then
                  LAPACK_LIBS=""
          fi
  fi
  
  # LAPACK linked to by default?  (is sometimes included in BLAS lib)
  if test $acx_lapack_ok = no; then
          save_LIBS="$LIBS"; LIBS="$LIBS $BLAS_LIBS $FLIBS"
          AC_CHECK_FUNC($cheev, [acx_lapack_ok=yes])
          LIBS="$save_LIBS"
  fi
  
  # Generic LAPACK library?
  for lapack in lapack lapack_rs6k; do
          if test $acx_lapack_ok = no; then
                  save_LIBS="$LIBS"; LIBS="$BLAS_LIBS $LIBS"
                  AC_CHECK_LIB($lapack, $cheev,
                      [acx_lapack_ok=yes; LAPACK_LIBS="-l$lapack"], [], [$FLIBS])
                  LIBS="$save_LIBS"
          fi
  done
  echo "LAPACK: " $LAPACK_LIBS
  echo "BLAS: " $BLAS_LIBS
  echo "LIBS: " $LIBS
  echo "FLIBS: " $FLIBS
  
  AC_SUBST(LAPACK_LIBS)
  
  # Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
  if test x"$acx_lapack_ok" = xyes; then
          ifelse([$1],,AC_DEFINE(HAVE_LAPACK,1,[Define if you have LAPACK library.]),[$1])
          :
  else
          acx_lapack_ok=no
          $2
  fi
  ])dnl ACX_LAPACK
