
include_directories(${Surfpack_SOURCE_DIR}/src)

## nkm library sources.
set(nkm_src
  NKM_GradKrigingModel.cpp
  NKM_KrigingModel.cpp
  NKM_LinearRegressionModel.cpp
  NKM_Optimize.cpp
  NKM_SurfData.cpp
  NKM_SurfPack.cpp
)

add_library(nkm ${nkm_src})
#target_link_libraries(surfaces ${LAPACK_LIBRARIES})
install(TARGETS nkm EXPORT ${Dakota_export_target} DESTINATION lib)

