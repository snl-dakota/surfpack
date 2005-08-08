dnl Checks whether GCC compiler is being used
AC_DEFUN([ACX_USING_GCC], [
  AC_PREREQ(2.50)
  AC_REQUIRE([AC_PROG_CXX])
  AM_CONDITIONAL([USING_GCC_COMPILER], [test "$ac_cv_cxx_compiler_gnu" = "yes"])
  ])dnl ACX_USING_GCC
