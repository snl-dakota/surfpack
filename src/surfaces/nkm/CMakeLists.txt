## nkm directory sources.
set(nkm_sources
  NKM_KrigingModel.cpp
  NKM_Optimize.cpp
  NKM_SurfData.cpp
  NKM_SurfPack.cpp
  NKM_SurfPackModel.cpp
  NKM_SurfMat.cpp
  NKM_PivotChol.f
)

# Optionally build the ValidateMain binary
if (Surfpack_NKM_Tests)
  set(SURFPACK_LIBS "surfpack_interpreter;surfpack;surfpack_fortran")
  set(SURFPACK_TPL_LIBS "ncsuopt;conmin")
  set(SURFPACK_SYSTEM_LIBS "${Boost_LIBRARIES};${LAPACK_LIBS};${BLAS_LIBS}")

  link_directories(${Surfpack_BINARY_DIR}/src/interpreter)
  add_executable(ValidateMain NKM_ValidateMain.cpp)
  target_link_libraries(ValidateMain ${SURFPACK_LIBS} ${SURFPACK_TPL_LIBS} 
    ${SURFPACK_SYSTEM_LIBS}
    )
endif()
