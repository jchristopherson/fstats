@PACKAGE_INIT@

include(CMakeFindDependencyMacro)
find_dependency(OpenMP COMPONENTS Fortran)

if (NOT BUILD_SHARED_LIBS)
  find_dependency(ferror QUIET)
  find_dependency(linalg QUIET)
  find_dependency(collections QUIET)
endif()

if(NOT TARGET "@PROJECT_NAME@::@PROJECT_NAME@")
  include("${CMAKE_CURRENT_LIST_DIR}/@PROJECT_NAME@-targets.cmake")
endif()