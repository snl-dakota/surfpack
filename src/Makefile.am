SUBDIRS = conmin kriging mars
lib_LTLIBRARIES = libsurfpack.la
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
			ann/client.h \
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
libsurfpack_la_LIBADD =	$(top_builddir)/src/mars/libmars.la \
			$(top_builddir)/src/kriging/libkriging.la \
			$(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)
INCLUDES = -I$(top_srcdir) -I$(top_srcdir)/src -I$(top_srcdir)/src/surfaces -I$(top_srcdir)/src/ann -I$(top_srcdir)/src/interpreter
CLEANFILES = *.bb *.bbg *.da *.gcov *.ii *.s
		
