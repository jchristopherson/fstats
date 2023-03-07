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