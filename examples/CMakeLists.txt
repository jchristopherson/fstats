# Get helper macros and functions
include("${PROJECT_SOURCE_DIR}/cmake/helper.cmake")

# Include Dependencies
add_subdirectory(fplot)
add_subdirectory(fortran_csv_module)

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
target_include_directories(allan_example PUBLIC ${fplot_INCLUDE_DIR})
target_link_libraries(allan_example ${fortran-csv-module_LIBRARY})
target_include_directories(allan_example PUBLIC ${fortran-csv-module_INCLUDE_DIR})

# Box-Muller Example
add_executable(box_muller_example box_muller_example.f90)
target_link_libraries(box_muller_example fstats ${fplot_LIBRARY})
target_include_directories(box_muller_example PUBLIC ${fplot_INCLUDE_DIR})

# Rejection Sampling Example
add_executable(rejection_sample_example rejection_sample_example.f90)
target_link_libraries(rejection_sample_example fstats ${fplot_LIBRARY})
target_include_directories(rejection_sample_example PUBLIC ${fplot_INCLUDE_DIR})

# LOWESS Example
add_executable(lowess_example lowess_example.f90)
target_link_libraries(lowess_example fstats ${fplot_LIBRARY})
target_include_directories(lowess_example PUBLIC ${fplot_INCLUDE_DIR})

# Bootstrap Example
add_executable(bootstrap_example bootstrap_example.f90)
target_link_libraries(bootstrap_example fstats ${fplot_LIBRARY})
target_include_directories(bootstrap_example PUBLIC ${fplot_INCLUDE_DIR})

# MCMC Regression Example
add_executable(mcmc_regression_example mcmc_regression_example.f90)
target_link_libraries(mcmc_regression_example fstats)
target_link_libraries(mcmc_regression_example ${fplot_LIBRARY})
target_include_directories(mcmc_regression_example PUBLIC ${fplot_INCLUDE_DIR})

# MCMC Regression Example 2
add_executable(mcmc_regression_example_2 mcmc_regression_example_2.f90)
target_link_libraries(mcmc_regression_example_2 fstats)
target_link_libraries(mcmc_regression_example_2 ${fplot_LIBRARY})
target_include_directories(mcmc_regression_example_2 PUBLIC ${fplot_INCLUDE_DIR})

# MCMC Regression Example 3
add_executable(mcmc_regression_example_3 mcmc_regression_example_3.f90)
target_link_libraries(mcmc_regression_example_3 fstats)
target_link_libraries(mcmc_regression_example_3 ${fplot_LIBRARY})
target_include_directories(mcmc_regression_example_3 PUBLIC ${fplot_INCLUDE_DIR})

# Distribution Example
add_executable(distribution_example distribution_example.f90)
target_link_libraries(distribution_example fstats)
target_link_libraries(distribution_example ${fplot_LIBRARY})
target_include_directories(distribution_example PUBLIC ${fplot_INCLUDE_DIR})

# Interpolation Example
add_executable(interpolation_example interpolation_example.f90)
target_link_libraries(interpolation_example fstats)
target_link_libraries(interpolation_example ${fplot_LIBRARY})
target_include_directories(interpolation_example PUBLIC ${fplot_INCLUDE_DIR})
