# Locate the source directory
set(dir ${CMAKE_CURRENT_SOURCE_DIR})

# Define the source files
set(FSTATS_SOURCES
    ${dir}/fstats.f90
)
set(FSTATS_SOURCES ${FSTATS_SOURCES} PARENT_SCOPE)