TESTS = srftest 
check_PROGRAMS = $(TESTS)
srftest_SOURCES = SurfPointTest.h \
		SurfPointTest.cpp \
		SurfDataTest.h \
		SurfDataTest.cpp \
		unittests.h \
		unittests.cpp \
		srftestmain.cpp 
srftest_CXXFLAGS = $(CPPUNIT_CFLAGS)
#srftest_LDFLAGS = -L$(HOME)/local/lib -L../src
srftest_LDADD = $(top_builddir)/src/libsurfpack.la 
INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/surfaces 
CLEANFILES = *.bb *.bbg *.da *.gcov *.ii *.s writePoint.sp writePoint.txt
