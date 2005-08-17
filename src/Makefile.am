if NEEDS_EXTRA_LAPACK_ROUTINES
extra_lapack_routines=$(top_srcdir)/src/kriging/dgecon.f \
		      $(top_srcdir)/src/kriging/dgelq2.f \
		      $(top_srcdir)/src/kriging/dgelqf.f \
		      $(top_srcdir)/src/kriging/dgels.f \
		      $(top_srcdir)/src/kriging/dgeqr2.f \
		      $(top_srcdir)/src/kriging/dgeqrf.f \
		      $(top_srcdir)/src/kriging/dgerq2.f \
		      $(top_srcdir)/src/kriging/dgerqf.f \
		      $(top_srcdir)/src/kriging/dgglse.f \
		      $(top_srcdir)/src/kriging/dggrqf.f \
		      $(top_srcdir)/src/kriging/dgetf2.f \
		      $(top_srcdir)/src/kriging/dgetrf.f \
		      $(top_srcdir)/src/kriging/dgetri.f \
		      $(top_srcdir)/src/kriging/dlabad.f \
		      $(top_srcdir)/src/kriging/dlacon.f \
		      $(top_srcdir)/src/kriging/dlamch.f \
		      $(top_srcdir)/src/kriging/dlange.f \
		      $(top_srcdir)/src/kriging/dlapy2.f \
		      $(top_srcdir)/src/kriging/dlarfb.f \
		      $(top_srcdir)/src/kriging/dlarf.f \
		      $(top_srcdir)/src/kriging/dlarfg.f \
		      $(top_srcdir)/src/kriging/dlarft.f \
		      $(top_srcdir)/src/kriging/dlascl.f \
		      $(top_srcdir)/src/kriging/dlaset.f \
		      $(top_srcdir)/src/kriging/dlassq.f \
		      $(top_srcdir)/src/kriging/dlaswp.f \
		      $(top_srcdir)/src/kriging/dlatrs.f \
		      $(top_srcdir)/src/kriging/dorm2r.f \
		      $(top_srcdir)/src/kriging/dorml2.f \
		      $(top_srcdir)/src/kriging/dormlq.f \
		      $(top_srcdir)/src/kriging/dormqr.f \
		      $(top_srcdir)/src/kriging/dormrq.f \
		      $(top_srcdir)/src/kriging/dormr2.f \
		      $(top_srcdir)/src/kriging/drscl.f \
		      $(top_srcdir)/src/kriging/dtrti2.f \
		      $(top_srcdir)/src/kriging/dtrtri.f \
		      $(top_srcdir)/src/kriging/ieeeck.f \
		      $(top_srcdir)/src/kriging/ilaenv.f \
		      $(top_srcdir)/src/kriging/lsame.f \
		      $(top_srcdir)/src/kriging/xerbla.f
else
extra_lapack_routines=
endif

AM_LFLAGS = -+ -o$(LEX_OUTPUT_ROOT).c
AM_YFLAGS = -d

if USING_GCC_COMPILER
BUILT_SOURCES = $(top_srcdir)/src/interpreter/lexer.cpp \
		$(top_srcdir)/src/interpreter/surfparse.cpp \
		$(top_srcdir)/src/interpreter/surfparse.h

standalone_only_sources=$(top_srcdir)/src/interpreter/lexer.lpp \
			$(top_srcdir)/src/interpreter/surfparse.ypp \
			$(top_srcdir)/src/interpreter/surfparse.h \
			$(top_srcdir)/src/interpreter/SurfpackParser.cpp \
			$(top_srcdir)/src/interpreter/SurfpackParser.h \
			$(top_srcdir)/src/interpreter/FlexWrapper.cpp \
			$(top_srcdir)/src/interpreter/FlexWrapper.h \
			$(top_srcdir)/src/interpreter/SurfpackInterpreter.cpp \
			$(top_srcdir)/src/interpreter/SurfpackInterpreter.h
else
standalone_only_sources=
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

libsurfpack_la_SOURCES=$(top_srcdir)/src/SurfPoint.cpp \
		       $(top_srcdir)/src/SurfPoint.h \
		       $(standalone_only_sources) \
		       $(top_srcdir)/src/interpreter/SurfpackParserArgs.cpp \
		       $(top_srcdir)/src/interpreter/SurfpackParserArgs.h \
		       $(top_srcdir)/src/SurfData.cpp \
		       $(top_srcdir)/src/SurfData.h \
		       $(top_srcdir)/src/Surface.cpp \
		       $(top_srcdir)/src/Surface.h \
		       $(top_srcdir)/src/SurfScaler.cpp \
		       $(top_srcdir)/src/SurfScaler.h \
		       $(top_srcdir)/src/AxesBounds.cpp \
		       $(top_srcdir)/src/AxesBounds.h \
		       $(top_srcdir)/src/SurfpackMatrix.tpl \
		       $(top_srcdir)/src/SurfpackMatrix.h \
		       $(top_srcdir)/src/surfpack_system_headers.h \
		       $(top_srcdir)/src/surfaces/KrigingSurface.cpp \
		       $(top_srcdir)/src/surfaces/KrigingSurface.h \
		       $(top_srcdir)/src/surfaces/KrigingCPPSurface.cpp \
		       $(top_srcdir)/src/surfaces/KrigingCPPSurface.h \
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
		       $(top_srcdir)/interface/SurfDataGenerator.cpp \
		       $(top_srcdir)/interface/SurfDataGenerator.h \
		       $(top_srcdir)/src/conmin/conmin.f \
		       $(top_srcdir)/src/mars/mars36_fort.f \
		       $(top_srcdir)/src/kriging/krig_model.f \
		       $(top_srcdir)/src/kriging/cpp_f77_conmin_cntl.f \
		       $(extra_lapack_routines)

## _LIBADD is automake's way of saying that these libraries must also
## be linked in whenever libsurfpack.la is linked against.  The
## top_builddir variable gets set by automake. LAPACK_LIBS, BLAS_LIBS,
## LIBS, and FLIBS are all set by the configure script.
libsurfpack_la_LIBADD = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)
if USING_GCC_COMPILER
  libsurfpack_la_LIBADD += $(LEXLIB)
endif

## Specifies additional places for the compiler to look for files that
## are requested via #include.
INCLUDES = -I$(top_srcdir) -I$(top_srcdir)/src -I$(top_srcdir)/src/surfaces \
	   -I$(top_srcdir)/src/ann -I$(top_srcdir)/src/interpreter

## The Makefiles produced by automake create a target named "clean"
## which removes, by default, the .o and .lo files that are normally
## created during a build.  The CLEANFILES automake variable allows
## the user to specify additional files that may need to be cleaned
## up.  Here, the .bb, .bbg, .da and .gcov files are extensions of
## files that may be created by the gcov profiler.  The *.ii and *.s
## files are temporary files that may be produced during C/C++
## compilation.
CLEANFILES = *.bb *.bbg *.da *.gcov *.ii *.s surfparse.h surfparse.cpp lexer.cpp
krig_time_test = krig_time
krig_time_SOURCES = $(top_srcdir)/system_tests/instances/time_krig_builds.cpp
krig_time_LDADD = $(top_builddir)/src/libsurfpack.la

bin_PROGRAMS = $(krig_time_test)

if USING_GCC_COMPILER
  bin_PROGRAMS += surfpack 
  surfpack_SOURCES = $(top_srcdir)/interface/main.cpp 
  surfpack_LDADD = $(top_builddir)/src/libsurfpack.la
  CLEANFILES += lexer.cpp surfparse.h surfparse.cpp $(BUILT_SOURCES)
endif
