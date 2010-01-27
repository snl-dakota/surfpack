## Process this file with automake to produce 'Makefile.in'.

## Initialize the SUBDIRS variable; it'll be added to below.
SUBDIRS = 

## Package subdirectories to recurse into.

# WJB: advise from Brian Adams: conmin is small so build it - 2 copies is fine
##if BUILD_CONMIN
SUBDIRS += conmin
##endif

## WJB - ToDo: add SVN external to teuchos
#if BUILD_TEUCHOS
#SUBDIRS += teuchos
#endif

## Conditionally configured subdirectories should _not_ be recursed into
## during 'make distclean' or 'make maintainer-clean'.
DIST_SUBDIRS = $(SUBDIRS)

## These files should be deleted on 'make maintainer-clean'.
MAINTAINERCLEANFILES = Makefile.in

