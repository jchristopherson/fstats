set(fstats_test_sources
    fstats_tests.f90
    fstats_test_helper.f90
    fstats_distribution_tests.f90
    fstats_statistics_tests.f90
    fstats_regression_tests.f90
    fstats_experimental_design_tests.f90
    fstats_nonlinear_regression_tests.f90
)

add_executable(fstats_tests ${fstats_test_sources})
target_link_libraries(fstats_tests fstats)
add_test(
    NAME fstats_tests
    WORKING_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
    COMMAND $<TARGET_FILE:fstats_tests>
)
