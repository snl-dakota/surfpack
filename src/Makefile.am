if NEEDS_EXTRA_LAPACK_ROUTINES
extra_lapack_routines=  $(top_srcdir)/src/kriging/dgecon.f \
			$(top_srcdir)/src/kriging/dlacon.f \
			$(top_srcdir)/src/kriging/dlamch.f \
			$(top_srcdir)/src/kriging/dlatrs.f \
			$(top_srcdir)/src/kriging/drscl.f
else
extra_lapack_routines=
endif

# lib_LTLIBRARIES is automake's way of saying "libraries that should be created
# and installed in the lib subdirectory of the top-level install directory. The
# "LT" in "LTLIBRARIES" means that libtool should be used to figure out whether
# a shared or static library is called for.
lib_LTLIBRARIES = libsurfpack.la

# libsurfpack_la is the canonicalized form of the libsurfpack.la that was 
# listed above as one of the targets that needs to be built by this makefile.
# (See the automake manual for the canonicalization rules.) The .h files are
# listed so that they will be included in the distribution, even though they
# may not need to be compiled. The _SOURCES part is automake's way of 
# specifying which files need to be compiled in order to build the 
# libsurfpack.la target.

libsurfpack_la_SOURCES =$(top_srcdir)/src/SurfPoint.cpp \
			$(top_srcdir)/src/SurfPoint.h \
			$(top_srcdir)/src/SurfData.cpp \
			$(top_srcdir)/src/SurfData.h \
			$(top_srcdir)/src/Surface.cpp \
			$(top_srcdir)/src/Surface.h \
			$(top_srcdir)/src/SurfScaler.cpp \
			$(top_srcdir)/src/SurfSclaer.h \
			$(top_srcdir)/src/AxesBounds.cpp \
			$(top_srcdir)/src/AxesBounds.h \
			$(top_srcdir)/src/surfaces/KrigingSurface.cpp \
			$(top_srcdir)/src/surfaces/KrigingSurface.h \
			$(top_srcdir)/src/surfaces/MarsSurface.cpp \
			$(top_srcdir)/src/surfaces/MarsSurface.h \
			$(top_srcdir)/src/surfaces/RBFNetSurface.cpp \
			$(top_srcdir)/src/surfaces/RBFNetSurface.h \
			$(top_srcdir)/src/surfaces/ANNSurface.cpp \
			$(top_srcdir)/src/surfaces/ANNSurface.h \
			$(top_srcdir)/src/surfaces/PolynomialSurface.cpp \
			$(top_srcdir)/src/surfaces/PolynomialSurface.h \
			$(top_srcdir)/src/surfpack.cpp \
			$(top_srcdir)/src/surfpack.h \
			$(top_srcdir)/src/SurfaceFactory.cpp \
			$(top_srcdir)/src/SurfaceFactory.h \
			$(top_srcdir)/src/ann/ann.cpp \
			$(top_srcdir)/src/ann/ann.h  \
			$(top_srcdir)/src/ann/client.cpp \
			$(top_srcdir)/src/ann/allocate.cpp \
			$(top_srcdir)/src/ann/allocate.h \
			$(top_srcdir)/src/ann/convert.cpp \
			$(top_srcdir)/src/ann/convert.h \
			$(top_srcdir)/src/ann/dsvd2.cpp \
			$(top_srcdir)/src/ann/dsvd2.h \
			$(top_srcdir)/src/ann/gen_drno.cpp \
			$(top_srcdir)/src/ann/gen_drno.h \
			$(top_srcdir)/src/ann/random.cpp \
			$(top_srcdir)/src/ann/random.h \
			$(top_srcdir)/src/ann/system_defs.h \
			$(top_srcdir)/src/ann/utilities.cpp \
			$(top_srcdir)/src/ann/utilities.h \
			$(top_srcdir)/src/ann/vector_enhancements.cpp \
			$(top_srcdir)/src/ann/vector_enhancements.h \
			$(top_srcdir)/src/interpreter/SurfpackParser.cpp \
			$(top_srcdir)/src/interpreter/SurfpackParser.h \
			$(top_srcdir)/src/interpreter/SurfpackInterpreter.cpp \
			$(top_srcdir)/src/interpreter/SurfpackInterpreter.h \
			$(top_srcdir)/src/interpreter/surfparse.tab.cpp \
			$(top_srcdir)/src/interpreter/surfparse.tab.h \
			$(top_srcdir)/src/interpreter/lex.yy.cpp \
			$(top_srcdir)/interface/SurfDataGenerator.cpp \
			$(top_srcdir)/interface/SurfDataGenerator.h \
			$(top_srcdir)/src/conmin/conmin.f \
			$(top_srcdir)/src/mars/mars36_fort.f \
			$(top_srcdir)/src/kriging/krig_model.f \
			$(top_srcdir)/src/kriging/cpp_f77_conmin_cntl.f \
			$(extra_lapack_routines)

# _LIBADD is automake's way of saying that these libraries must also be linked
# in whenever libsurfpack.la is linked against.  The top_builddir variable
# gets set by automake. LAPACK_LIBS, BLAS_LIBS, LIBS, and FLIBS are all set
# by the configure script.

libsurfpack_la_LIBADD = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)

# Specifies additional places for the compiler to look for files that are
# requested via #include.
INCLUDES = -I$(top_srcdir) -I$(top_srcdir)/src -I$(top_srcdir)/src/surfaces \
	   -I$(top_srcdir)/src/ann -I$(top_srcdir)/src/interpreter

# The Makefiles produces by automake create a target named "clean" which 
# removes, by default, the .o and .lo files that are normally created during
# a build.  The CLEANFILES automake variable allows the user to specify
# additional files that may need to be cleaned up.  Here, the .bb,  .bbg,
# .da and .gcov files are extensions of files that may be created by the
# gcov profiler.  The *.ii and *.s files are temporary fils that may be
# produce during C/C++ compilation. 
CLEANFILES = *.bb *.bbg *.da *.gcov *.ii *.s
		
bin_PROGRAMS = surfpack 
surfpack_SOURCES = $(top_srcdir)/interface/main.cpp 
surfpack_LDADD = $(top_builddir)/src/libsurfpack.la 
