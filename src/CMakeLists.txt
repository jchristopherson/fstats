# Locate the source directory
set(dir ${CMAKE_CURRENT_SOURCE_DIR})

# Define the source files
set(FSTATS_SOURCES
    ${dir}/fstats.f90
    ${dir}/fstats_special_functions.f90
    ${dir}/fstats_descriptive_statistics.f90
    ${dir}/fstats_hypothesis.f90
    ${dir}/fstats_anova.f90
    ${dir}/fstats_distributions.f90
    ${dir}/fstats_helper_routines.f90
    ${dir}/fstats_regression.f90
    ${dir}/fstats_experimental_design.f90
    ${dir}/fstats_allan.f90
    ${dir}/fstats_errors.f90
    ${dir}/fstats_bootstrap.f90
    ${dir}/fstats_sampling.f90
    ${dir}/fstats_smoothing.f90
    ${dir}/fstats_types.f90
    ${dir}/fstats_mcmc.f90
    ${dir}/fstats_interp.f90
)
set(FSTATS_SOURCES ${FSTATS_SOURCES} PARENT_SCOPE)

# Enable preprocessor directives
set_source_files_properties(
    ${FSTATS_SOURCES}
    PROPERTIES Fortran_PREPROCESS ON
)

# Dependencies
include(FetchContent)
find_package(OpenMP REQUIRED COMPONENTS Fortran)
find_package(ferror QUIET)
find_package(linalg QUIET)
find_package(collections QUIET)
if (NOT ferror_FOUND)
    message(STATUS "FERROR could not be found.  A reference version will be employed.")
    FetchContent_Declare(
        ferror
        GIT_REPOSITORY "https://github.com/jchristopherson/ferror"
    )
    FetchContent_MakeAvailable(ferror)
    set(FERROR_LIBRARIES ferror)
    set(BUILD_FERROR TRUE)
else()
    set(FERROR_LIBRARIES ferror::ferror)
    set(BUILD_FERROR FALSE)
endif()

if (NOT linalg_FOUND)
    message(STATUS "LINALG could not be found.  A reference version will be employed.")
    FetchContent_Declare(
        linalg
        GIT_REPOSITORY "https://github.com/jchristopherson/linalg"
    )
    FetchContent_MakeAvailable(linalg)
    set(LINALG_LIBRARIES linalg)
    set(BUILD_LINALG TRUE)
else()
    set(LINALG_LIBRARIES linalg::linalg)
    set(BUILD_LINALG FALSE)
endif()

if (NOT collections_FOUND)
    message(STATUS "COLLECTIONS could not be found.  A reference version will be employed.")
    FetchContent_Declare(
        collections
        GIT_REPOSITORY "https://github.com/jchristopherson/collections"
        GIT_TAG main
    )
    FetchContent_MakeAvailable(collections)
    set(COLLECTIONS_LIBRARIES collections)
    set(BUILD_COLLECTIONS TRUE)
else()
    set(COLLECTIONS_LIBRARIES collections::collections)
    set(BUILD_COLLECTIONS FALSE)
endif()

# Build
add_library(${PROJECT_NAME} ${FSTATS_SOURCES})
target_link_libraries(
    ${PROJECT_NAME}
    PUBLIC
    ${FERROR_LIBRARIES}
    ${COLLECTIONS_LIBRARIES}
)
target_link_libraries(
    ${PROJECT_NAME}
    PRIVATE
    ${LINALG_LIBRARIES}
)
if (OpenMP_Fortran_FOUND)
    target_link_libraries(${PROJECT_NAME} PRIVATE OpenMP::OpenMP_Fortran)
    target_compile_definitions(${PROJECT_NAME} PRIVATE USEOPENMP=1)
endif()

set_target_properties(
    ${PROJECT_NAME}
    PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    WINDOWS_EXPORT_ALL_SYMBOLS ON
)

set(FSTATS_MOD_DIR ${CMAKE_CURRENT_BINARY_DIR}/mod_files/)
if (NOT EXISTS "${FSTATS_MOD_DIR}")
    file(MAKE_DIRECTORY "${FSTATS_MOD_DIR}")
endif()

set_target_properties(
    ${PROJECT_NAME} PROPERTIES
    Fortran_MODULE_DIRECTORY ${FSTATS_MOD_DIR}
)
target_include_directories(
    ${PROJECT_NAME} PUBLIC
    $<BUILD_INTERFACE:${FSTATS_MOD_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_MODULEDIR}>
)

# Install
install(
    TARGETS ${PROJECT_NAME}
    EXPORT ${PROJECT_NAME}-targets
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
install(DIRECTORY ${FSTATS_MOD_DIR} DESTINATION ${CMAKE_INSTALL_MODULEDIR})

if (${BUILD_FERROR} AND ${BUILD_SHARED_LIBS})
    install(IMPORTED_RUNTIME_ARTIFACTS ${FERROR_LIBRARIES})
endif()

if (${BUILD_LINALG} AND ${BUILD_SHARED_LIBS})
    install(IMPORTED_RUNTIME_ARTIFACTS ${LINALG_LIBRARIES})
endif()

if (${BUILD_COLLECTIONS} AND ${BUILD_SHARED_LIBS})
    install(IMPORTED_RUNTIME_ARTIFACTS ${COLLECTIONS_LIBRARIES})
endif()
