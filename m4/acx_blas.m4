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
