# Automake if/else/endif construction.  The MDR_CXX_CPPUNIT macro that is
# invoked in autoconf invokes and AM_CONDITIONAL macro for INCLUDE_TESTS.
# In other words, if the test to link against CPP Unit succeeds, the
# "if INCLUDE_TESTS" will be true. 

# If CPP Unit was found, descend into the tests subdirectory and build
# the unit tests.
if INCLUDE_TESTS
test_dirs=tests
endif

# Specifies that make should also descend into these subdirectories during the
# build process.
SUBDIRS = src $(test_dirs)

# Specifies that the files in the examples subdirectory should be included in 
# the distribution tarball, even though they are not involved in the build
# process.
EXTRA_DIST = $(top_srcdir)/examples/* $(top_srcdir)/docs/DoxygenConfig
