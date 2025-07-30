cmake_minimum_required(VERSION 3.26)

project(fstats Fortran)
enable_language(Fortran)

set(CMAKE_Fortran_STANDARD 2018)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)

set(MT_LIBS_DIR "${CMAKE_BINARY_DIR}/../../../libs")

add_subdirectory(configure)
add_subdirectory(dependencies)

find_package(OpenMP REQUIRED COMPONENTS Fortran)

set(FSTATS_SOURCES
    src/fstats.f90
    src/fstats_special_functions.f90
    src/fstats_descriptive_statistics.f90
    src/fstats_hypothesis.f90
    src/fstats_anova.f90
    src/fstats_distributions.f90
    src/fstats_helper_routines.f90
    src/fstats_regression.f90
    src/fstats_experimental_design.f90
    src/fstats_allan.f90
    src/fstats_errors.f90
    src/fstats_bootstrap.f90
    src/fstats_sampling.f90
    src/fstats_smoothing.f90
    src/fstats_types.f90
    src/fstats_mcmc.f90
    src/fstats_interp.f90
)
set_source_files_properties(
    ${FSTATS_SOURCES}
    PROPERTIES
    Fortran_PREPROCESS ON
)

set(MT_SOURCES
    mt/mt_fstats.f90
)

add_library(fstats SHARED
    ${FSTATS_SOURCES}
    ${MT_SOURCES}
)

set_target_properties(
    fstats
    PROPERTIES
    OUTPUT_NAME fstats
    POSITION_INDEPENDENT_CODE TRUE
    LIBRARY_OUTPUT_DIRECTORY ${MT_LIBS_DIR}
)
target_include_directories(
    fstats
    PUBLIC
    $<BUILD_INTERFACE:${mod_dir}>
    $<BUILD_INTERFACE:${ferror_INCLUDE_DIR}>
    $<BUILD_INTERFACE:${linalg_INCLUDE_DIR}>
    $<BUILD_INTERFACE:${collections_INCLUDE_DIR}>
    $<INSTALL_INTERFACE:${include_install_dir}>
)

target_link_directories(fstats PUBLIC
    ${MT_LIBS_DIR}
)

target_link_libraries(fstats PUBLIC
    ${ferror_LIBRARY} 
    ${linalg_LIBRARY} 
    ${collections_LIBRARY}
)

if (OpenMP_Fortran_FOUND)
    target_link_libraries(fstats PUBLIC OpenMP::OpenMP_Fortran)
    target_compile_definitions(fstats PRIVATE USEOPENMP=1)
endif()

target_compile_options(fstats PRIVATE
    $<$<CONFIG:DEBUG>:-O0 -g>
    $<$<CONFIG:RELEASE>:-O3 -s>
    -std=f2018
    -Wall -Werror -pedantic -pedantic-errors
    -Wno-error=maybe-uninitialized
    -Wno-error=unused-variable
    -Wno-error=unused-dummy-argument
    -Wno-error=unused-label
    -Wno-error=surprising
    -Wno-error=uninitialized
    -fimplicit-none
)

add_custom_command(TARGET fstats POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
    $<TARGET_FILE:fstats>
    ${MT_LIBS_DIR}
)
if (WIN32)
    add_custom_command(TARGET fstats POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_BINARY_DIR}/_deps/lapack-build/bin/libblas.dll
        ${CMAKE_BINARY_DIR}/_deps/lapack-build/bin/liblapack.dll
        ${MT_LIBS_DIR}
    )
endif()
