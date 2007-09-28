## Process this file with automake to produce Makefile.in

AM_YFLAGS = -d

if BUILD_STANDALONE
noinst_LTLIBRARIES = libinterpreter.la
libinterpreter_la_SOURCES = lexer.l \
			     surfparse.ypp \
			     surfparse.h \
			     SurfpackParser.cpp \
			     SurfpackParser.h \
			     FlexWrapper.cpp \
			     FlexWrapper.h \
			     SurfpackInterpreter.cpp \
			     SurfpackInterpreter.h

BUILT_SOURCES = lexer.c \
		surfparse.cpp \
		surfparse.h
endif

INCLUDES = -I$(top_srcdir) -I$(top_srcdir)/src

MAINTAINERCLEANFILES = Makefile.in

CLEANFILES = $(BUILT_SOURCES) *.bb *.bbg *.da *.gcov *.ii *.s
