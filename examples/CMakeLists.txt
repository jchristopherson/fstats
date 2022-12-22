# Full Factorial Example
add_executable(full_factorial_example full_factorial_example.f90)
target_link_libraries(full_factorial_example fstats)

# ANOVA Example
add_executable(anova_example_1 anova_example_1.f90)
target_link_libraries(anova_example_1 fstats)