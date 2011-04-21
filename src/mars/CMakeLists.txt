## Mars library sources.
set(mars_SOURCES
  mars36_fort.f		
)

## Set the CONMIN library name.
add_library(mars ${mars_SOURCES})
#set_target_properties(mars PROPERTIES LINKER_LANGUAGE Fortran)
install(TARGETS mars EXPORT ${ExportTarget} DESTINATION lib)

