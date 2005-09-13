## Process this file with automake to produce Makefile.in
noinst_LTLIBRARIES = libinterpreter.la

libinterpreter_la_SOURCES = SurfpackParserArgs.cpp \
			    SurfpackParserArgs.h

AM_LFLAGS = -o$(LEX_OUTPUT_ROOT).c
AM_YFLAGS = -d

if BUILD_STANDALONE
libinterpreter_la_SOURCES += lexer.l \
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

CLEANFILES = *.bb *.bbg *.da *.gcov *.ii *.s
MAINTAINERCLEANFILES = Makefile.in
