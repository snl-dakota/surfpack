bin_PROGRAMS = surfpack
surfpack_SOURCES = main.cpp
surfpack_LDADD = $(top_srcdir)/src/libsurfpack.la \
		 $(top_srcdir)/src/ann/libann.la \
		 $(top_srcdir)/src/mars/libmars.la \
		 $(top_srcdir)/src/kriging/libkriging.la \
		 $(top_srcdir)/src/conmin/libconmin.la \
		 -llapack -lblas -lm -lg2c
INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/surfaces 
