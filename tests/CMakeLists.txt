# Get helper macros and functions
include("${PROJECT_SOURCE_DIR}/cmake/helper.cmake")

# Dependencies
add_subdirectory(fortran_test_helper)

# Define the source files
set(fstats_test_sources
    fstats_tests.f90
    fstats_test_helper.f90
    fstats_distribution_tests.f90
    fstats_statistics_tests.f90
    fstats_regression_tests.f90
    fstats_experimental_design_tests.f90
    fstats_nonlinear_regression_tests.f90
    fstats_allan_tests.f90
    fstats_mcmc_tests.f90
)

# Build the tests
add_executable(fstats_tests ${fstats_test_sources})
link_library(fstats_tests ${PROJECT_NAME} ${PROJECT_INCLUDE_DIR})
link_library(fstats_tests ${fortran_test_helper_LIBRARY} ${fortran_test_helper_INCLUDE_DIR})
add_test(
    NAME fstats_tests
    WORKING_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
    COMMAND $<TARGET_FILE:fstats_tests>
)
