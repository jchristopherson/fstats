# Locate the source directory
set(dir ${CMAKE_CURRENT_SOURCE_DIR})

# Define the source files
set(FSTATS_SOURCES
    ${dir}/fstats.f90
    ${dir}/fstats_special_functions.f90
    ${dir}/fstats_descriptive_statistics.f90
    ${dir}/fstats_hypothesis.f90
    ${dir}/fstats_anova.f90
    ${dir}/fstats_distributions.f90
    ${dir}/fstats_helper_routines.f90
    ${dir}/fstats_regression.f90
    ${dir}/fstats_experimental_design.f90
    ${dir}/fstats_allan.f90
    ${dir}/fstats_errors.f90
    ${dir}/fstats_bootstrap.f90
    ${dir}/fstats_sampling.f90
    ${dir}/fstats_smoothing.f90
    ${dir}/fstats_types.f90
    ${dir}/fstats_mcmc.f90
    # ${dir}/fstats_mcmc_fitting.f90
)
set(FSTATS_SOURCES ${FSTATS_SOURCES} PARENT_SCOPE)
