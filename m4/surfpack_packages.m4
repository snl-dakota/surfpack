dnl Packages

AC_DEFUN([SURFPACK_PACKAGES],[
  AC_ARG_WITH([conmin],
              AC_HELP_STRING([--with-conmin=DIR],
                             [use CONMIN library in specified DIR]),
              [],[with_conmin="yes"])
  dnl acx_local_conmin=no
  acx_local_conmin=yes
  case $with_conmin in
  no)
    dnl CONMIN package is needed unconditionally.
    dnl Surfpack provides source code for CONMIN in the packages directory.
    AC_MSG_ERROR([Surfpack cannot be configured without CONMIN. Please specify
                  directory to conmin OR simply, --with-conmin=yes
                  to get the default path to the Surfpack provided conmin pkg.])
    ;;

  dnl For yes, check CONMIN_ROOT, otherwise fallback to local CONMIN
  yes | "")
    AC_MSG_CHECKING([for Conmin])
    if test -n "$CONMIN_ROOT" -a -d "$CONMIN_ROOT"; then

      AC_MSG_RESULT([using Conmin in CONMIN_ROOT: $CONMIN_ROOT])

    elif test -d `pwd`/packages/conmin; then

      dnl Use local Conmin and instruct subpackages to do so as well
      export CONMIN_ROOT=`pwd`/packages/conmin
      acx_local_conmin=yes

      dnl Configure the local conmin package for use in Surfpack
      dnl WJB: temporarily, build conmin in all cases -- AC_CONFIG_SUBDIRS([packages/conmin])
      AC_MSG_RESULT([using local Conmin in $CONMIN_ROOT])

    else
      AC_MSG_NOTICE([could not find Conmin directory.])
      AC_MSG_NOTICE([need help locating conmin!])
      AC_MSG_ERROR([PLEASE PROVIDE full path to conmin, --with-conmin=<DIR>])
    fi
    ;;

  dnl Otherwise, user should have provided an explicit path to Conmin
  *)
    AC_MSG_CHECKING([for specified Conmin])
    CONMIN_ROOT=$withval
    if test -n "$CONMIN_ROOT" -a -d "$CONMIN_ROOT"; then
      AC_MSG_RESULT([using: $CONMIN_ROOT])
    else
      AC_MSG_ERROR([could not locate $CONMIN_ROOT])
    fi
    ;;

  esac

  dnl WJB: temporarily, build conmin in all cases
  AC_CONFIG_SUBDIRS([packages/conmin])
  CONMIN_LDFLAGS="-L$CONMIN_ROOT/src"
  AC_ARG_VAR(CONMIN_ROOT, [Path to Conmin, an Optimization library in written in F77])
  AC_SUBST(CONMIN_LDFLAGS)

  AM_CONDITIONAL([BUILD_CONMIN], [test "x$acx_local_conmin" = xyes])


  dnl -------------------------
  dnl Teuchos include DIR check
  dnl Teuchos lib DIR check
  dnl -------------------------
  dnl WJB - ToDo: copy/paste from DAKOTA's packages M4 file
  dnl Surfpack depends on Teuchos UNCONDITIONALLY (not quite yet)
dnl  AC_ARG_WITH([teuchos-lib],
dnl              AC_HELP_STRING([--with-teuchos-lib=DIR],
dnl                             [use Teuchos libraries in specified lib DIR]),
dnl              [],[with_teuchos_lib="yes"])

dnl  AC_ARG_VAR(TEUCHOS_ROOT, [Path to Teuchos, OO Numerics foundation library])
dnl  AC_SUBST(TEUCHOS_CPPFLAGS)
dnl  AC_SUBST(TEUCHOS_LDFLAGS)

dnl  AM_CONDITIONAL([BUILD_TEUCHOS], [test "x$acx_local_teuchos" = xyes])
])

