bin_PROGRAMS = surfpack
surfpack_SOURCES = main.cpp
surfpack_LDADD = $(top_builddir)/src/libsurfpack.la \
		 $(top_builddir)/src/ann/libann.la \
		 $(top_builddir)/src/mars/libmars.la \
		 $(top_builddir)/src/kriging/libkriging.la \
		 $(top_builddir)/src/conmin/libconmin.la \
		 -llapack -lblas -lg2c
INCLUDES = -I$(top_srcdir) -I$(top_srcdir)/src -I$(top_srcdir)/src/surfaces 
CLEANFILES = *.bb *.bbg *.da *.gcov *.ii *.s
