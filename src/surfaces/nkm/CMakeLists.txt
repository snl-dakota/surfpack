
include_directories(${PROJECT_SOURCE_DIR}/src)

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
#TARGET_LINK_LIBRARIES(surfaces ${LAPACK_LIBRARIES})
INSTALL(TARGETS nkm DESTINATION lib)

