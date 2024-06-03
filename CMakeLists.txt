# Master CMAKE Build Script
cmake_minimum_required(VERSION 3.24)
project(
    fstats
    LANGUAGES Fortran
    VERSION 1.2.4
)

# Confgiure everything
add_subdirectory(configure)

# Deal with the dependencies
add_subdirectory(dependencies)

# Get OpenMP
find_package(OpenMP REQUIRED COMPONENTS Fortran)

# Source
add_subdirectory(src)
set_source_files_properties(
    ${FSTATS_SOURCES}
    PROPERTIES Fortran_PREPROCESS ON
)
add_fortran_library(
    ${PROJECT_NAME}
    ${PROJECT_INCLUDE_DIR}
    ${CMAKE_INSTALL_INCLUDEDIR}
    ${PROJECT_VERSION}
    ${PROJECT_VERSION_MAJOR}
    ${FSTATS_SOURCES}
)
target_link_libraries(${PROJECT_NAME} PUBLIC ${ferror_LIBRARY} ${linalg_LIBRARY})
target_include_directories(${PROJECT_NAME} PUBLIC
    $<BUILD_INTERFACE:${ferror_INCLUDE_DIR}>
    $<BUILD_INTERFACE:${linalg_INCLUDE_DIR}>)
if (OpenMP_Fortran_FOUND)
    target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_Fortran)
    target_compile_definitions(${PROJECT_NAME} PRIVATE USEOPENMP=1)
endif()

# Testing
option(BUILD_TESTING "Build tests")
include(CTest)
message(STATUS "Build tests: ${BUILD_TESTING}")
if (BUILD_TESTING)
    enable_testing()
    add_subdirectory(tests)
endif()

# Examples
option(BUILD_FSTATS_EXAMPLES "Build FSTATS examples")
message(STATUS "Build FSTATS examples: ${BUILD_FSTATS_EXAMPLES}")
if (BUILD_FSTATS_EXAMPLES)
    add_subdirectory(examples)
endif()