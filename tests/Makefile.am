TESTS = srftest 
check_PROGRAMS = $(TESTS)
srftest_SOURCES = SurfPointTest.h SurfPointTest.cpp SurfDataTest.h SurfDataTest.cpp srftestmain.cpp 
srftest_CXXFLAGS = $(CPPUNIT_CFLAGS)
srftest_LDFLAGS = -L$(HOME)/local/lib -L../src
srftest_LDADD = $(top_builddir)/src/libsurfpack.la \
		 $(top_builddir)/src/ann/libann.la \
		 $(top_builddir)/src/mars/libmars.la \
		 $(top_builddir)/src/kriging/libkriging.la \
		 $(top_builddir)/src/conmin/libconmin.la \
		 -lcppunit -ldl \
		 -llapack -lblas -lm -lg2c
INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/surfaces 
CLEANFILES = *.bb *.bbg *.da *.gcov
