lib_LTLIBRARIES = libkriging.la
if NEEDS_EXTRA_LAPACK_ROUTINES
extra_lapack_routines= dgecon.f dlacon.f dlamch.f dlatrs.f drscl.f
else
extra_lapack_routines=
endif

libkriging_la_SOURCES = krig_model.f \
			cpp_f77_conmin_cntl.f \
			$(extra_lapack_routines)

libkriging_la_LIBADD = $(top_builddir)/src/conmin/libconmin.la $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)
CLEANFILES = *.bb *.bbg *.da *.gcov *.ii *.s
