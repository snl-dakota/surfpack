bin_PROGRAMS = surfpack 
#polybench datawriter driv
surfpack_SOURCES = main.cpp \
		   SurfDataGenerator.cpp \
		   SurfDataGenerator.h
surfpack_LDADD = $(top_builddir)/src/libsurfpack.la 
#polybench_SOURCES = polybench.cpp $(top_srcdir)/tests/unittests.h $(top_srcdir)/tests/unittests.cpp
#datawriter_SOURCES = datawriter.cpp SurfDataGenerator.cpp SurfDataGenerator.h
#driv_SOURCES = drivmain.cpp 
#polybench_LDADD = $(top_builddir)/src/libsurfpack.la
#datawriter_LDADD = $(top_builddir)/src/libsurfpack.la 
#driv_LDADD = $(top_builddir)/src/libsurfpack.la 

INCLUDES = -I$(top_srcdir) -I$(top_srcdir)/src -I$(top_srcdir)/src/surfaces -I$(top_srcdir)/tests -I$(top_srcdir)/src/interpreter
CLEANFILES = *.bb *.bbg *.da *.gcov *.ii *.s
