dnl Checks whether a suitable LEX and YACC have been found
dnl If so, sets a flag to create the Surfpack standalone 
AC_DEFUN([ACX_BUILD_STANDALONE], [
  AC_PREREQ(2.50)
  AC_REQUIRE([AC_PROG_LEX])
  AC_REQUIRE([AC_PROG_YACC])
  AM_CONDITIONAL([BUILD_STANDALONE], [test "$ac_cv_prog_LEX" != "" -a "$ac_cv_prog_YACC" != ""])
  ])dnl ACX_BUILD_STANDALONE
