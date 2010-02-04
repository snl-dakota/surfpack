
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/src)

## surfaces library sources.
SET(surfaces_src
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

ADD_LIBRARY(surfaces ${surfaces_src})
#TARGET_LINK_LIBRARIES(surfaces ${LAPACK_LIBRARIES})

