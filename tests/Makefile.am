# Specifies a list of executables that should be built and run when the user
# types "make check."  If TESTS is not defined and the names of the executables
# are listed directly under check_PROGRAMS.  Those executables will be built
# on a "make check" command but will not execute automatically.
TESTS = srftest
check_PROGRAMS = $(TESTS) 

# These are the source files needed for the srftest target.
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

# Somehow, CPP Unit is able to communicate its desired C++ compilation flags.
# Frankly, I don't know how this works.
srftest_CXXFLAGS = $(CPPUNIT_CFLAGS)

# Specifies flags for the srftest link line.
srftest_LDFLAGS = -L$(top_builddir)/../../lib 

# Specifies what libraries are necessary for the srftest build
# LAPACK_LIBS, BLAS_LIBS, and FLIBS would be set by the configure script.
srftest_LDADD = $(top_builddir)/src/libsurfpack.la $(CPPUNIT_LIBS) \
		$(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)

# Specifies additional places to look for files requested via #include
INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/surfaces \
	   -I$(top_srcdir)/src/interpreter


# The Makefiles produces by automake create a target named "clean" which 
# removes, by default, the .o and .lo files that are normally created during
# a build.  The CLEANFILES automake variable allows the user to specify
# additional files that may need to be cleaned up.  Here, the .bb,  .bbg,
# .da and .gcov files are extensions of files that may be created by the
# gcov profiler.  The *.ii and *.s files are temporary fils that may be
# produce during C/C++ compilation. 
CLEANFILES = *.bb *.bbg *.da *.gcov *.ii *.s
