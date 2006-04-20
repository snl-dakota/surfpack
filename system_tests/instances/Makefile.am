## Process this file with automake to produce Makefile.in
krig_time_test = krig_time
krig_time_SOURCES = time_krig_builds.cpp
krig_time_LDADD = $(top_builddir)/src/libsurfpack.la
if NEEDS_EXTRA_LAPACK_ROUTINES
krig_time_LDADD += $(BLAS_LIBS) $(LIBS) $(FLIBS) $(LEXLIB)
else
krig_time_LDADD += $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS) $(LEXLIB)
endif

bin_PROGRAMS = $(krig_time_test)

INCLUDES = -I$(top_srcdir) -I$(top_srcdir)/src \
	   -I$(top_srcdir)/src/interpreter \
	   -I$(top_srcdir)/src/surfaces

CLEANFILES = *.bb *.bbg *.da *.gcov *.ii *.s

MAINTAINERCLEANFILES = Makefile.in
