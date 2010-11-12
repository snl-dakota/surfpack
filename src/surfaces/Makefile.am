## Process this file with automake to produce Makefile.in
SUBDIRS = nkm

SUBLIBS = nkm/libnkm.la

noinst_LTLIBRARIES = libsurfaces.la

libsurfaces_la_SOURCES =  MovingLeastSquaresModel.cpp \
			 MovingLeastSquarsModel.h \
			 LinearRegressionModel.cpp \
			 LinearRegressionModel.h \
			 RadialBasisFunctionModel.cpp \
			 RadialBasisFunctionModel.h \
			 DirectANNModel.cpp \
			 DirectANNModel.h \
		         KrigingModel.cpp \
			 KrigingModel.h \
		         MarsModel.cpp \
			 MarsModel.h \
			 SurfpackModel.cpp \
			 SurfpackModel.h \
			 ModelScaler.cpp \
			 ModelScaler.h

libsurfaces_la_LIBADD = $(SUBLIBS)

INCLUDES = -I$(top_srcdir) -I$(top_srcdir)/src \
	   -I$(top_srcdir)/src/interpreter \
	   -I$(top_srcdir)/src/surfaces

CLEANFILES = *.bb *.bbg *.da *.gcov *.ii *.s
MAINTAINERCLEANFILES = Makefile.in
