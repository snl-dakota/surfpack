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
