bin_PROGRAMS = surfpack polybench
surfpack_SOURCES = main.cpp
polybench_SOURCES = polybench.cpp $(top_srcdir)/tests/unittests.h $(top_srcdir)/tests/unittests.cpp
polybench_LDADD = $(top_builddir)/src/libsurfpack.la
surfpack_LDADD = $(top_builddir)/src/libsurfpack.la 
INCLUDES = -I$(top_srcdir) -I$(top_srcdir)/src -I$(top_srcdir)/src/surfaces -I$(top_srcdir)/tests
CLEANFILES = *.bb *.bbg *.da *.gcov *.ii *.s
