# Locate the source directory
set(dir ${CMAKE_CURRENT_SOURCE_DIR})

# Define the source files
set(FSTATS_SOURCES
    ${dir}/fstats.f90
    ${dir}/special_functions_beta.f90
    ${dir}/special_functions_gamma.f90
    ${dir}/distributions_normal.f90
    ${dir}/distributions_t.f90
    ${dir}/distributions_f.f90
    ${dir}/statistics_implementation.f90
    ${dir}/regression_implementation.f90
    ${dir}/experimental_design_implementation.f90
    ${dir}/levenberg_marquardt.f90
    ${dir}/distributions_chisquared.f90
    ${dir}/allan.f90
    ${dir}/fstats_errors.f90
    ${dir}/external.f90
    ${dir}/bootstrapping.f90
    ${dir}/distributions_binomial.f90
)
set(FSTATS_SOURCES ${FSTATS_SOURCES} PARENT_SCOPE)