include_directories(${Surfpack_SOURCE_DIR}/src)

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
#target_link_libraries(surfaces ${LAPACK_LIBRARIES})
install(TARGETS surfaces DESTINATION lib)

