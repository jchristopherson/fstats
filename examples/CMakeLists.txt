# Include Dependencies
add_subdirectory(fortran_csv_module)
include(FetchContent)
FetchContent_Declare(
    fplot
    GIT_REPOSITORY "https://github.com/jchristopherson/fplot"
)
FetchContent_MakeAvailable(fplot)
set(fplot_LIBRARY fplot)

# Full Factorial Example
add_executable(full_factorial_example full_factorial_example.f90)
target_link_libraries(full_factorial_example fstats)

# ANOVA Example
add_executable(anova_example_1 anova_example_1.f90)
target_link_libraries(anova_example_1 fstats)

# Regression Example
add_executable(regression_example regression_example.f90)
target_link_libraries(regression_example fstats)

# Nonlinear Regression Example
add_executable(nl_regression_example nl_regression_example.f90)
target_link_libraries(nl_regression_example fstats)

# Allan Variance Example
configure_file(
    ${CMAKE_SOURCE_DIR}/examples/data/noise_data.csv
    ${CMAKE_BINARY_DIR}/examples/noise_data.csv
    COPYONLY
)
add_executable(allan_example allan_example.f90)
target_link_libraries(allan_example fstats)
target_link_libraries(allan_example ${fplot_LIBRARY})
target_link_libraries(allan_example ${fortran-csv-module_LIBRARY})
target_include_directories(allan_example PUBLIC ${fortran-csv-module_INCLUDE_DIR})

# Box-Muller Example
add_executable(box_muller_example box_muller_example.f90)
target_link_libraries(box_muller_example fstats ${fplot_LIBRARY})

# Rejection Sampling Example
add_executable(rejection_sample_example rejection_sample_example.f90)
target_link_libraries(rejection_sample_example fstats ${fplot_LIBRARY})

# LOWESS Example
add_executable(lowess_example lowess_example.f90)
target_link_libraries(lowess_example fstats ${fplot_LIBRARY})

# Bootstrap Example
add_executable(bootstrap_example bootstrap_example.f90)
target_link_libraries(bootstrap_example fstats ${fplot_LIBRARY})

# MCMC Regression Example
add_executable(mcmc_regression_example mcmc_regression_example.f90)
target_link_libraries(mcmc_regression_example fstats)
target_link_libraries(mcmc_regression_example ${fplot_LIBRARY})

# MCMC Regression Example 2
add_executable(mcmc_regression_example_2 mcmc_regression_example_2.f90)
target_link_libraries(mcmc_regression_example_2 fstats)
target_link_libraries(mcmc_regression_example_2 ${fplot_LIBRARY})

# MCMC Regression Example 3
add_executable(mcmc_regression_example_3 mcmc_regression_example_3.f90)
target_link_libraries(mcmc_regression_example_3 fstats)
target_link_libraries(mcmc_regression_example_3 ${fplot_LIBRARY})

# Distribution Example
add_executable(distribution_example distribution_example.f90)
target_link_libraries(distribution_example fstats)
target_link_libraries(distribution_example ${fplot_LIBRARY})

# Interpolation Example
add_executable(interpolation_example interpolation_example.f90)
target_link_libraries(interpolation_example fstats)
target_link_libraries(interpolation_example ${fplot_LIBRARY})
