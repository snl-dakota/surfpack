TESTS = srftest 
check_PROGRAMS = $(TESTS)
srftest_SOURCES = SurfPointTest.h SurfPointTest.cpp srftestmain.cpp 
srftest_CXXFLAGS = $(CPPUNIT_CFLAGS)
srftest_LDFLAGS = -L$(HOME)/local/lib -L../src
srftest_LDADD = $(top_srcdir)/src/libsurfpack.la \
		 $(top_srcdir)/src/ann/libann.la \
		 $(top_srcdir)/src/mars/libmars.la \
		 $(top_srcdir)/src/kriging/libkriging.la \
		 $(top_srcdir)/src/conmin/libconmin.la \
		 -lcppunit -ldl \
		 -llapack -lblas -lm -lg2c
INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/surfaces 
