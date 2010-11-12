## Process this file with automake to produce Makefile.in
noinst_LTLIBRARIES = libnkm.la

libnkm_la_SOURCES =  NKM_KrigingModel.cpp \
	NKM_KrigingModel.hpp \
	NKM_LinearRegressionModel.cpp \
	NKM_LinearRegressionModel.hpp \
	NKM_Optimize.cpp \
	NKM_Optimize.hpp \
	NKM_SurfData.cpp \
	NKM_SurfData.hpp \
	NKM_SurfMat.hpp \
	NKM_SurfMat_Native.hpp \
	NKM_SurfPack.cpp \
	NKM_SurfPack.hpp \
	NKM_SurfPackModel.hpp
#	NKM_SurfMat_Teuchos.hpp

bin_PROGRAMS = TestMain ValidateMain

TestMain_SOURCES = NKM_TestMain.cpp
TestMain_LDADD = libnkm.la
TestMain_LDADD += $(top_builddir)/packages/CONMIN/src/libconmin.la
TestMain_LDADD += $(top_builddir)/packages/NCSUOpt/libncsuopt.a
TestMain_LDADD += $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)

ValidateMain_SOURCES = NKM_ValidateMain.cpp
ValidateMain_LDADD = libnkm.la
ValidateMain_LDADD += $(top_builddir)/packages/CONMIN/src/libconmin.la
ValidateMain_LDADD += $(top_builddir)/packages/NCSUOpt/libncsuopt.a
ValidateMain_LDADD += $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)

MAINTAINERCLEANFILES = Makefile.in
