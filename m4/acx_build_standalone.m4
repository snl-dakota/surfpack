dnl Checks whether a suitable LEX and YACC have been found
dnl If so, sets a flag to create the Surfpack standalone 
AC_DEFUN([ACX_BUILD_STANDALONE], [
  AC_PREREQ(2.50)
  AC_REQUIRE([AC_PROG_LEX])
  AC_REQUIRE([AC_PROG_YACC])
  AC_CACHE_CHECK([whether to build standalone executable 'surfpack'],
    [ax_cv_build_standalone],
    [if test -n "$ac_cv_prog_LEX"; dnl add additional logic to check for yacc
       then ax_cv_build_standalone=yes; 
       else ax_cv_build_standalone=no; 
     fi]) 
  AM_CONDITIONAL([BUILD_STANDALONE], [test $ax_cv_build_standalone=yes])
  ])dnl ACX_BUILD_STANDALONE
