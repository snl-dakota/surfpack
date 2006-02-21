## Process this file with automake to produce Makefile.in

## Add all the necessary m4 directories to aclocal's search path.
ACLOCAL_AMFLAGS = -I m4

## Specifies that make should also descend into these subdirectories
## during the build process.
SUBDIRS = src interface system_tests tests

## Specifies that the files in the examples subdirectory should be
## included in the distribution tarball, even though they are not
## involved in the build process.
EXTRA_DIST = $(top_srcdir)/examples/* $(top_srcdir)/docs/DoxygenConfig

## The Makefiles produced by automake create a target named "clean"
## which removes, by default, the .o and .lo files that are normally
## created during a build.  The CLEANFILES automake variable allows the
## user to specify additional files that may need to be cleaned up.
## Here, the .bb, .bbg, .da and .gcov files are extensions of files
## that may be created by the gcov profiler.  The *.ii and *.s files
## are temporary files that may be produced during C/C++ compilation.
CLEANFILES = *.bb *.bbg *.da *.gcov *.ii *.s

## MAINTAINERCLEANFILES serves a similar purpose to CLEANFILES, but
## is only used for the maintainer-clean target.
MAINTAINERCLEANFILES = Makefile.in src/Makefile.in tests/Makefile.in \
		       aclocal.m4 compile config.guess config.sub \
		       configure depcomp install-sh ltmain.sh missing \
		       surfpack_config.h.in tests/badscale.txt
