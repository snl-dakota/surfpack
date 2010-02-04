## Mars library sources.
SET(mars_SOURCES
  mars36_fort.f		
)

## Set the CONMIN library name.
ADD_LIBRARY(mars ${mars_SOURCES})
#set_target_properties(mars PROPERTIES LINKER_LANGUAGE Fortran)

