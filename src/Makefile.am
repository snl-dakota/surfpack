# Specifies subdirectories into which the make utility should descend during
# the build process.
SUBDIRS = conmin kriging mars

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

libsurfpack_la_SOURCES = SurfPoint.cpp \
			SurfPoint.h \
			SurfData.cpp \
			SurfData.h \
			Surface.cpp \
			Surface.h \
			PointDefinition.cpp \
			PointDefinition.h \
			surfaces/KrigingSurface.cpp \
			surfaces/KrigingSurface.h \
			surfaces/MarsSurface.cpp \
			surfaces/MarsSurface.h \
			surfaces/RBFNetSurface.cpp \
			surfaces/RBFNetSurface.h \
			surfaces/ANNSurface.cpp \
			surfaces/ANNSurface.h \
			surfaces/PolynomialSurface.cpp \
			surfaces/PolynomialSurface.h \
			surfpack.cpp \
			surfpack.h \
			SurfaceFactory.cpp \
			SurfaceFactory.h \
			ann/ann.cpp \
			ann/ann.h  \
			ann/client.cpp \
			ann/allocate.cpp \
			ann/allocate.h \
			ann/convert.cpp \
			ann/convert.h \
			ann/dsvd2.cpp \
			ann/dsvd2.h \
			ann/gen_drno.cpp \
			ann/gen_drno.h \
			ann/random.cpp \
			ann/random.h \
			ann/system_defs.h \
			ann/utilities.cpp \
			ann/utilities.h \
			ann/vector_enhancements.cpp \
			ann/vector_enhancements.h \
			interpreter/SurfpackParser.cpp \
			interpreter/SurfpackParser.h \
			interpreter/SurfpackInterpreter.cpp \
			interpreter/SurfpackInterpreter.h \
			interpreter/surfparse.tab.cpp \
			interpreter/surfparse.tab.h \
			interpreter/lex.yy.cpp

# _LIBADD is automake's way of saying that these libraries must also be linked
# in whenever libsurfpack.la is linked against.  The top_builddir variable
# gets set by automake. LAPACK_LIBS, BLAS_LIBS, LIBS, and FLIBS are all set
# by the configure script.

libsurfpack_la_LIBADD =	$(top_builddir)/src/mars/libmars.la \
			$(top_builddir)/src/kriging/libkriging.la \
			$(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)

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
		
