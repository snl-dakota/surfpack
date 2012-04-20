include_directories("${Surfpack_SOURCE_DIR}/src" 
  "${Surfpack_SOURCE_DIR}/src/interpreter"
)

set(SURFPACK_LIBS "surfpack_interpreter;surfpack")
set(SURFPACK_TPL_LIBS "ncsuopt;conmin")
set(SURFPACK_SYSTEM_LIBS "${FLEX_LIBRARIES};${LAPACK_LIBS}")

link_directories(${Surfpack_BINARY_DIR}/src/interpreter)
add_executable(surfpack_exe main.cpp)
target_link_libraries(surfpack_exe ${SURFPACK_LIBS} ${SURFPACK_TPL_LIBS} 
  ${SURFPACK_SYSTEM_LIBS}
)
