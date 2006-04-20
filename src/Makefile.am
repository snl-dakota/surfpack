## Process this file with automake to produce Makefile.in
SUBDIRS = ann conmin interpreter kriging mars surfaces

SUBLIBS = ann/libann.la \
          conmin/libconmin.la \
          interpreter/libinterpreter.la \
          kriging/libkriging.la \
          mars/libmars.la \
          surfaces/libsurfaces.la

if NEEDS_EXTRA_LAPACK_ROUTINES
SUBDIRS += lapack
SUBLIBS += lapack/liblapack.la
endif

## lib_LTLIBRARIES is automake's way of saying "libraries that should
## be created and installed in the lib subdirectory of the top-level
## install directory. The "LT" in "LTLIBRARIES" means that libtool
## should be used to figure out whether a shared or static library is
## called for.
lib_LTLIBRARIES = libsurfpack.la

## libsurfpack_la is the canonicalized form of the libsurfpack.la that
## was listed above as one of the targets that needs to be built by
## this makefile.  (See the automake manual for the canonicalization
## rules.) The .h files are listed so that they will be included in
## the distribution, even though they may not need to be compiled. The
## _SOURCES part is automake's way of specifying which files need to
## be compiled in order to build the libsurfpack.la target.
libsurfpack_la_SOURCES = SurfPoint.cpp \
			 SurfPoint.h \
			 SurfData.cpp \
			 SurfData.h \
			 Surface.cpp \
			 Surface.h \
			 SurfScaler.cpp \
			 SurfScaler.h \
			 AxesBounds.cpp \
			 AxesBounds.h \
			 SurfpackMatrix.tpl \
			 SurfpackMatrix.h \
			 surfpack_system_headers.h \
			 surfpack.cpp \
			 surfpack.h \
			 SurfaceFactory.cpp \
			 SurfaceFactory.h

## _LIBADD is automake's way of saying that these libraries must also
## be linked in whenever libsurfpack.la is linked against.  The
## top_builddir variable gets set by automake.
libsurfpack_la_LIBADD = $(SUBLIBS)

## Specifies additional places for the compiler to look for files that
## are requested via #include.
INCLUDES = -I$(top_srcdir) -I$(top_srcdir)/src \
	   -I$(top_srcdir)/src/ann \
	   -I$(top_srcdir)/src/interpreter \
	   -I$(top_srcdir)/src/surfaces

## The Makefiles produced by automake create a target named "clean"
## which removes, by default, the .o and .lo files that are normally
## created during a build.  The CLEANFILES automake variable allows
## the user to specify additional files that may need to be cleaned
## up.  Here, the .bb, .bbg, .da and .gcov files are extensions of
## files that may be created by the gcov profiler.  The *.ii and *.s
## files are temporary files that may be produced during C/C++
## compilation.
CLEANFILES = *.bb *.bbg *.da *.gcov *.ii *.s
