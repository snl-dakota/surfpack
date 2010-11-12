dnl Packages

AC_DEFUN([SURFPACK_PACKAGES],[

  AC_CONFIG_SUBDIRS([packages/CONMIN])
  AC_CONFIG_SUBDIRS([packages/NCSUOpt])

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

