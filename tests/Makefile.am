TESTS = srftest 
check_PROGRAMS = $(TESTS)
srftest_SOURCES = SurfPointTest.h \
		SurfPointTest.cpp \
		unittests.h \
		unittests.cpp \
		SurfDataTest.h \
		SurfDataTest.cpp \
		PolynomialSurfaceTest.cpp \
		PolynomialSurfaceTest.h \
		SurfaceTest.cpp \
		SurfaceTest.h \
		srftestmain.cpp 
srftest_CXXFLAGS = $(CPPUNIT_CFLAGS)
srftest_LDFLAGS = -L$(top_builddir)/../../lib 
srftest_LDADD = $(top_builddir)/src/libsurfpack.la $(CPPUNIT_LIBS) $(top_builddir)/src/mars/libmars.la $(top_builddir)/src/kriging/libkriging.la $(top_builddir)/src/conmin/libconmin.la $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)
INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/surfaces -I$(top_srcdir)/src/interpreter
CLEANFILES = *.bb *.bbg *.da *.gcov *.ii *.s
