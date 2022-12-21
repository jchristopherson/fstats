# Master CMAKE Build Script
cmake_minimum_required(VERSION 3.17)
project(
    fstats
    LANGUAGES Fortran C
    VERSION 1.0.0
)

# Utilize the GNU directory structure
include(GNUInstallDirs)

# Get helper macros and functions
include("${PROJECT_SOURCE_DIR}/cmake/helper.cmake")

# Confgiure everything
add_subdirectory(configure)

# Deal with the dependencies
find_package(ferror 1.4.0 QUIET)
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
link_library(${PROJECT_NAME} ${ferror_LIBRARY} ${ferror_INCLUDE_DIR})

# Installation

# Testing
option(BUILD_TESTING "Build tests")
include(CTest)
message(STATUS "Build tests: ${BUILD_TESTING}")
if (BUILD_TESTING)
    enable_testing()
endif()

# Examples
option(BUILD_FSTATS_EXAMPLES "Build FSTATS examples")
message(STATUS "Build FSTATS examples: ${BUILD_FSTATS_EXAMPLES}")
if (BUILD_FSTATS_EXAMPLES)
endif()