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

# Bootstrap Regression Example
add_executable(bootstrap_regression_example bootstrap_regression_example.f90)
target_link_libraries(bootstrap_regression_example fstats)