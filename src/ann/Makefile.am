## Process this file with automake to produce Makefile.in
noinst_LTLIBRARIES = libann.la

libann_la_SOURCES = ann.cpp \
                    ann.h \
                    allocate.cpp \
                    allocate.h \
                    convert.cpp \
                    convert.h \
                    dsvd2.cpp \
                    dsvd2.h \
                    gen_drno.cpp \
                    gen_drno.h \
                    random.cpp \
                    random.h \
                    utilities.cpp \
                    utilities.h \
                    vector_enhancements.cpp \
                    vector_enhancements.h

INCLUDES = -I$(top_srcdir) -I$(top_srcdir)/src -I$(top_srcdir)/src/ann

CLEANFILES = *.bb *.bbg *.da *.gcov *.ii *.s
MAINTAINERCLEANFILES = Makefile.in
