@PACKAGE_INIT@

include(CMakeFindDependencyMacro)
find_dependency(OpenMP COMPONENTS Fortran)
find_dependency(linalg 3.0 QUIET)
find_dependency(collections 1.2.1 QUIET)
find_dependency(fstring 1.1 QUIET)

if(NOT TARGET "@PROJECT_NAME@::@PROJECT_NAME@")
  include("${CMAKE_CURRENT_LIST_DIR}/@PROJECT_NAME@-targets.cmake")
endif()