# TODO: consider checking in the parser files and possibly configuring
# into the build tree, conditionally using source vs. binary.

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

bison_target(SurfParser surfparse.ypp ${CMAKE_CURRENT_BINARY_DIR}/surfparse.cpp)
flex_target(SurfScanner lexer.l  ${CMAKE_CURRENT_BINARY_DIR}/lexer.c)
add_flex_bison_dependency(SurfScanner SurfParser)

set(surfparse_hpp "${CMAKE_CURRENT_BINARY_DIR}/surfparse.hpp")
set(surfparse_h "${CMAKE_CURRENT_BINARY_DIR}/surfparse.h")

add_custom_command(
  DEPENDS "${surfparse_hpp}"
  OUTPUT  "${surfparse_h}"
  COMMAND ${CMAKE_COMMAND} -E copy "${surfparse_hpp}" "${surfparse_h}"
)

set(surfpack_interpreter_sources
  ${CMAKE_CURRENT_BINARY_DIR}/lexer.c
  ${CMAKE_CURRENT_BINARY_DIR}/surfparse.cpp
  ${CMAKE_CURRENT_BINARY_DIR}/surfparse.h
  SurfpackParser.cpp
  SurfpackParser.h
  FlexWrapper.cpp
  FlexWrapper.h
  SurfpackInterpreter.cpp
  SurfpackInterpreter.h
)

add_library(surfpack_interpreter ${surfpack_interpreter_sources})
