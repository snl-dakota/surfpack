## Process this file with automake to produce Makefile.in

## Build the Surfpack binary iff BUILD_STANDALONE has been set 
## Must be explicitly enabled at configure time
if BUILD_STANDALONE
  bin_PROGRAMS = surfpack
  surfpack_SOURCES = main.cpp 
  surfpack_LDADD = $(top_builddir)/src/libsurfpack.la
  surfpack_LDADD += $(top_builddir)/packages/CONMIN/src/libconmin.la
  surfpack_LDADD += $(top_builddir)/packages/NCSUOpt/libncsuopt.a
if HAVE_BOOST
  surfpack_CPPFLAGS = @BOOST_CPPFLAGS@
  surfpack_LDFLAGS = @BOOST_LDFLAGS@
  surfpack_LDADD += -lboost_serialization
endif
  # Shouldn't need ${LEXLIB} any more due to %noyywrap in lexer.l
  surfpack_LDADD += $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)
endif


INCLUDES = -I$(top_srcdir) -I$(top_srcdir)/src \
	   -I$(top_srcdir)/src/interpreter

CLEANFILES = *.bb *.bbg *.da *.gcov *.ii *.s

MAINTAINERCLEANFILES = Makefile.in
