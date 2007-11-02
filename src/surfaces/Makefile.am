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
			 PolynomialSurface.h \
			 MovingLeastSquaresModel.cpp \
			 MovingLeastSquarsModel.h \
			 LinearRegressionModel.cpp \
			 LinearRegressionModel.h \
			 RadialBasisFunctionModel.cpp \
			 RadialBasisFunctionModel.h \
			 DirectANNModel.cpp \
			 DirectANNModel.h \
		         KrigingModel.cpp \
			 KrigingModel.h \
			 SurfpackModel.cpp \
			 SurfpackModel.h \
			 ModelScaler.cpp \
			 ModelScaler.h

INCLUDES = -I$(top_srcdir) -I$(top_srcdir)/src \
	   -I$(top_srcdir)/src/ann \
	   -I$(top_srcdir)/src/interpreter \
	   -I$(top_srcdir)/src/surfaces

CLEANFILES = *.bb *.bbg *.da *.gcov *.ii *.s
MAINTAINERCLEANFILES = Makefile.in
