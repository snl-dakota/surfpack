SUBDIRS = ann conmin kriging mars 
lib_LTLIBRARIES = libsurfpack.la
libsurfpack_la_SOURCES = SurfPoint.cpp \
			SurfPoint.h \
			SurfData.cpp \
			SurfData.h \
			Surface.cpp \
			Surface.h \
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
			surfpack.h
INCLUDES = -I$(top_srcdir)/src/surfaces -I$(top_srcdir)/src/ann
		
