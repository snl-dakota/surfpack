## Process this file with automake to produce Makefile.in
##-- SUBDIRS = ann interpreter kriging mars surfaces
SUBDIRS = mars surfaces
#SUBDIRS = ann kriging mars surfaces

SUBLIBS = mars/libmars.la \
          surfaces/libsurfaces.la

if BUILD_STANDALONE
SUBDIRS += interpreter 
SUBLIBS += interpreter/libinterpreter.la
endif

## lib_LTLIBRARIES is automake's way of saying "libraries that should
## be created and installed in the lib subdirectory of the top-level
## install directory. The "LT" in "LTLIBRARIES" means that libtool
## should be used to figure out whether a shared or static library is
## called for.
lib_LTLIBRARIES = libsurfpack.la

##libsurfpack_la_CPPFLAGS= -g -fno-inline
libsurfpack_la_CPPFLAGS=
## libsurfpack_la is the canonicalized form of the libsurfpack.la that
## was listed above as one of the targets that needs to be built by
## this makefile.  (See the automake manual for the canonicalization
## rules.) The .h files are listed so that they will be included in
## the distribution, even though they may not need to be compiled. The
## _SOURCES part is automake's way of specifying which files need to
## be compiled in order to build the libsurfpack.la target.

libsurfpack_la_SOURCES = AxesBounds.cpp \
			 AxesBounds.h \
			 MersenneTwister.h \
			 ModelFactory.cpp \
			 ModelFactory.h \
			 SurfData.cpp \
			 SurfData.h \
			 surfpack.cpp \
			 surfpack.h \
			 SurfpackMatrix.h \
			 SurfpackMatrix.tpl \
			 surfpack_system_headers.h \
			 SurfpackInterface.h \
			 SurfpackInterface.cpp \
			 SurfpackParserArgs.h \
                         SurfpackParserArgs.cpp \
			 SurfPoint.cpp \
			 SurfPoint.h \
			 Conmin.cpp \
			 Comin.h \
			 ModelFitness.cpp \
			 ModelFitness.h

## Given the near-term move to CMake, only distributing headers we know 
## external clients depend on
include_HEADERS = SurfpackMatrix.h \
	          surfpack_system_headers.h

## _LIBADD is automake's way of saying that these libraries must also
## be linked in whenever libsurfpack.la is linked against.  The
## top_builddir variable gets set by automake.
libsurfpack_la_LIBADD = $(SUBLIBS)

## Specifies additional places for the compiler to look for files that
## are requested via #include.
INCLUDES = -I$(top_builddir) -I$(top_srcdir) -I$(top_srcdir)/src \
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
