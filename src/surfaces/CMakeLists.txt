
include_directories(${PROJECT_SOURCE_DIR}/src)

subdirs(nkm)

## surfaces library sources.
set(surfaces_src
   MovingLeastSquaresModel.cpp
   MovingLeastSquaresModel.h
   LinearRegressionModel.cpp
   LinearRegressionModel.h
   RadialBasisFunctionModel.cpp
   RadialBasisFunctionModel.h
   DirectANNModel.cpp
   DirectANNModel.h
   KrigingModel.cpp
   KrigingModel.h
   MarsModel.cpp
   MarsModel.h
   SurfpackModel.cpp
   SurfpackModel.h
   ModelScaler.cpp
   ModelScaler.h
)

add_library(surfaces ${surfaces_src})
#TARGET_LINK_LIBRARIES(surfaces ${LAPACK_LIBRARIES})
INSTALL(TARGETS surfaces DESTINATION lib)

