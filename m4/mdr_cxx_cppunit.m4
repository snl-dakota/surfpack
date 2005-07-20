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
