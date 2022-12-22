# Master CMAKE Build Script
cmake_minimum_required(VERSION 3.24)
project(
    fstats
    LANGUAGES Fortran
    VERSION 1.0.0
)

# Utilize the GNU directory structure
include(GNUInstallDirs)

# Get helper macros and functions
include("${PROJECT_SOURCE_DIR}/cmake/helper.cmake")

# Confgiure everything
add_subdirectory(configure)

# Deal with the dependencies
find_package(BLAS)
find_package(LAPACK)
find_package(ferror 1.4.0 QUIET)
find_package(linalg 1.7.0 QUIET)
add_subdirectory(dependencies)

# Source
add_subdirectory(src)
add_fortran_library(
    ${PROJECT_NAME}
    ${PROJECT_INCLUDE_DIR}
    ${CMAKE_INSTALL_INCLUDEDIR}
    ${PROJECT_VERSION}
    ${PROJECT_VERSION_MAJOR}
    ${FSTATS_SOURCES}
)
target_link_libraries(${PROJECT_NAME} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
link_library(${PROJECT_NAME} ${ferror_LIBRARY} ${ferror_INCLUDE_DIR})
link_library(${PROJECT_NAME} ${linalg_LIBRARY} ${linalg_INCLUDE_DIR})

# Installation
add_subdirectory(install)

# Testing
option(BUILD_FSTATS_TESTING "Build tests")
include(CTest)
message(STATUS "Build tests: ${BUILD_FSTATS_TESTING}")
if (BUILD_FSTATS_TESTING)
    enable_testing()
    add_subdirectory(tests)
endif()

# Examples
option(BUILD_FSTATS_EXAMPLES "Build FSTATS examples")
message(STATUS "Build FSTATS examples: ${BUILD_FSTATS_EXAMPLES}")
if (BUILD_FSTATS_EXAMPLES)
    add_subdirectory(examples)
endif()