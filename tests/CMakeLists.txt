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
include(FetchContent)
find_package(fortran_test_helper QUIET)
find_package(ferror QUIET)
find_package(linalg QUIET)
if (NOT fortran_test_helper_FOUND)
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
if (NOT ferror_FOUND)
    FetchContent_Declare(
        ferror
        GIT_TAG "origin/master"
        GIT_REPOSITORY https://github.com/jchristopherson/ferror
    )
    FetchContent_MakeAvailable(ferror)
    set(FERROR_LIBRARIES ferror)
else()
    set(FERROR_LIBRARIES ferror::ferror)
endif()
if (NOT linalg_FOUND)
    FetchContent_Declare(
        linalg
        GIT_TAG "origin/master"
        GIT_REPOSITORY https://github.com/jchristopherson/linalg
    )
    FetchContent_MakeAvailable(linalg)
    set(LINALG_LIBRARIES linalg)
else()
    set(LINALG_LIBRARIES linalg::linalg)
endif()

# Build the tests
add_executable(fstats_tests ${fstats_test_sources})
target_link_libraries(
    fstats_tests
    fstats
    ${FTH_LIBRARIES}
    ${FERROR_LIBRARIES}
    ${LINALG_LIBRARIES}
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

