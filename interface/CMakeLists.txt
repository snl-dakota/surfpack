include_directories("${Surfpack_SOURCE_DIR}/src" 
  "${Surfpack_SOURCE_DIR}/src/interpreter"
)

set(SURFPACK_LIBS "surfpack_interpreter;surfpack;surfpack_fortran")
set(SURFPACK_TPL_LIBS "ncsuopt;conmin")
# Shouldn't need ${FLEX_LIBRARIES} any more due to %noyywrap in lexer.l
set(SURFPACK_SYSTEM_LIBS "${Boost_LIBRARIES};${LAPACK_LIBS};${BLAS_LIBS}")


link_directories(${Surfpack_BINARY_DIR}/src/interpreter)
add_executable(surfpack_exe main.cpp)
target_link_libraries(surfpack_exe ${SURFPACK_LIBS} ${SURFPACK_TPL_LIBS} 
  ${SURFPACK_SYSTEM_LIBS}
)
set_target_properties(surfpack_exe PROPERTIES OUTPUT_NAME surfpack) 

install(TARGETS surfpack_exe DESTINATION bin)
