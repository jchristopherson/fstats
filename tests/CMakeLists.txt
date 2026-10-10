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
    fstats_bootstrap_tests.f90
    fstats_missing_data_tests.f90
    fstats_msa_tests.f90
    fstats_glm_tests.f90
)

# Dependencies
include(FetchContent)
find_package(fortran_test_helper QUIET)
if (NOT fortran_test_helper_FOUND)
    set(BUILD_TESTING OFF)
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

add_test(NAME fstats_glm COMMAND $<TARGET_FILE:fstats_tests> glm)
foreach(invalid_case IN ITEMS
    rows columns empty no-columns underdetermined family nan-x inf-y nan-beta
    binomial-range poisson-range gamma-range iterations tolerance nan-tolerance
    cutoff inf-cutoff singular zero-column gamma-domain poisson-overflow
    poisson-underflow binomial-boundary predictor-overflow sum-overflow
    negative nan infinite size unallocated zero stats-size covariance-size
    alpha-zero alpha-one alpha-nan gamma-inference-dof robust-inference-dof covariance-overflow
)
    add_test(NAME fstats_glm_invalid_${invalid_case}
        COMMAND $<TARGET_FILE:fstats_tests> glm-invalid-${invalid_case})
    set_tests_properties(fstats_glm_invalid_${invalid_case} PROPERTIES WILL_FAIL TRUE)
endforeach()

foreach(thread_count IN ITEMS 1 4)
    add_test(
        NAME fstats_mcmc_${thread_count}_threads
        COMMAND $<TARGET_FILE:fstats_tests> mcmc
    )
    set_tests_properties(fstats_mcmc_${thread_count}_threads PROPERTIES
        ENVIRONMENT "OMP_NUM_THREADS=${thread_count};OMP_MAX_ACTIVE_LEVELS=2"
    )
endforeach()

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

