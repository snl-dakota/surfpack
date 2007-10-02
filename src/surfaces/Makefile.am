## Process this file with automake to produce Makefile.in
noinst_LTLIBRARIES = libsurfaces.la

libsurfaces_la_SOURCES = KrigingSurface.cpp \
			 KrigingSurface.h \
			 MarsSurface.cpp \
			 MarsSurface.h \
			 RBFNetSurface.cpp \
			 RBFNetSurface.h \
			 ANNSurface.cpp \
			 ANNSurface.h \
			 PolynomialSurface.cpp \
			 PolynomialSurface.h

INCLUDES = -I$(top_srcdir) -I$(top_srcdir)/src \
	   -I$(top_srcdir)/src/ann \
	   -I$(top_srcdir)/src/interpreter \
	   -I$(top_srcdir)/src/surfaces

CLEANFILES = *.bb *.bbg *.da *.gcov *.ii *.s
MAINTAINERCLEANFILES = Makefile.in
