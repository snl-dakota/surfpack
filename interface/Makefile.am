## Process this file with automake to produce Makefile.in

## Build the Surfpack binary iff BUILD_STANDALONE has been set 
if BUILD_STANDALONE
  bin_PROGRAMS = surfpack
  surfpack_SOURCES = main.cpp \
		     SurfDataGenerator.cpp \
		     SurfDataGenerator.h
  surfpack_LDADD = $(top_builddir)/src/libsurfpack.la
endif

INCLUDES = -I$(top_srcdir) -I$(top_srcdir)/src \
	   -I$(top_srcdir)/src/ann \
	   -I$(top_srcdir)/src/interpreter

CLEANFILES = *.bb *.bbg *.da *.gcov *.ii *.s
MAINTAINERCLEANFILES = Makefile.in
