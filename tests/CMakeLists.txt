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
    fstats_interp_tests.f90
)

# Dependencies
find_package(fortran_test_helper QUIET)
if (NOT fortran_test_helper_FOUND)
    include(FetchContent)
    FetchContent_Declare(
        fortran_test_helper
        GIT_TAG "origin/main"
        GIT_REPOSITORY https://github.com/jchristopherson/fortran_test_helper
    )
    FetchContent_MakeAvailable(fortran_test_helper)
    set(FTH_LIBRARIES fortran_test_helper)
else()
    set(FTH_LIBRARIES fortran_test_helper::fortran_test_helper)
endif()

# Build the tests
add_executable(fstats_tests ${fstats_test_sources})
target_link_libraries(
    fstats_tests
    fstats
    ${FTH_LIBRARIES}
)
add_test(
    NAME fstats_tests
    WORKING_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
    COMMAND $<TARGET_FILE:fstats_tests>
)

if (${BUILD_SHARED_LIBS} AND WIN32)
    add_custom_command(
        TARGET fstats_tests
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy 
        $<TARGET_RUNTIME_DLLS:fstats_tests> 
        $<TARGET_FILE_DIR:fstats_tests>
        COMMAND_EXPAND_LISTS
    )
endif()

