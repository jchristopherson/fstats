# fstats
FSTATS is a modern Fortran 2018 statistical library. The public API is collected in the `fstats` module, so applications can generally start with `use fstats`.

## Status
[![CMake](https://github.com/jchristopherson/fstats/actions/workflows/cmake.yml/badge.svg)](https://github.com/jchristopherson/fstats/actions/workflows/cmake.yml)
[![Actions Status](https://github.com/jchristopherson/fstats/workflows/fpm/badge.svg)](https://github.com/jchristopherson/fstats/actions)

## Capabilities
FSTATS includes the following areas of functionality:

- Descriptive statistics: means, variance, standard deviation, medians, quantiles, trimmed means, covariance, and pooled variance.
- Probability distributions: normal, log-normal, Student's t, F, chi-squared, binomial, Poisson, and multivariate normal distributions.
- Hypothesis testing: confidence intervals, t-tests, F-tests, Bartlett's test, Levene's test, and sample-size calculations.
- ANOVA: one-factor and two-factor analysis of variance.
- Regression: polynomial linear least squares, regression statistics, R-squared metrics, correlations, numerical Jacobians, and nonlinear Levenberg-Marquardt least squares.
- Experimental design: full and fractional factorial designs, central composite designs, Latin hypercube designs, model fitting, diagnostics, prediction, model comparison, ANOVA, efficiency, and response-surface optimization.
- Measurement systems analysis: gauge repeatability and reproducibility (gauge R&R) studies of crossed, nested, and expanded designs, reporting variance components, percent contribution, percent study variation, percent tolerance, and the number of distinct categories.
- Resampling and simulation: bootstrap resampling, random sampling, rejection sampling, Box-Muller sampling, and multivariate normal sampling.
- Markov chain Monte Carlo: chains, target distributions, proposals, samplers, and model evaluation.
- Signal and numerical methods: Allan variance, LOWESS smoothing, linear/polynomial/spline/Hermite interpolation, and missing-data imputation.
- Special functions: beta and gamma functions, regularized and incomplete forms, and the digamma function.

## Building with CMake
[CMake](https://cmake.org/) 3.24 or newer, a Fortran 2018 compiler, Git, and OpenMP support are required. CMake looks for [LINALG](https://github.com/jchristopherson/linalg) and [COLLECTIONS](https://github.com/jchristopherson/collections); if compatible installations are not found, it fetches reference versions automatically.

Configure and build the library:

```text
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
```

Tests and examples are disabled by default. Enable either or both at configure time:

```text
cmake -S . -B build -DBUILD_TESTING=ON -DBUILD_FSTATS_EXAMPLES=ON
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

Install the library, module files, CMake package files, and pkg-config metadata with:

```text
cmake --install build --prefix ./install
```

[FPM](https://github.com/fortran-lang/fpm) can also be used to build this library using the provided fpm.toml.
```text
fpm build --profile release
fpm test --profile release
```

FPM resolves the `linalg`, `collections`, and OpenMP dependencies declared by this project. To use FSTATS as a dependency in another FPM project, add the following to that project's `fpm.toml`:
```toml
[dependencies]
fstats = { git = "https://github.com/jchristopherson/fstats", tag = "v1.7.0" }
```

Then import the API in Fortran with `use fstats` and build normally with `fpm build`.

## Parallel MCMC
`sample_chains` runs independent Metropolis-Hastings chains with OpenMP. Supply
equally sized arrays of samplers, proposals, and targets. Each target must be
initialized independently, rather than copied from an object whose collections
or pointer components might share mutable state. The observations are shared
read-only; models and sampler callbacks must also be thread-safe.

For a user-defined `custom_target` extending `mcmc_target`, with `xdata` and
`ydata` already defined, declare and initialize the chain objects as follows:

```fortran
type(mcmc_sampler) :: samplers(4)
type(mcmc_proposal) :: proposals(4)
type(custom_target) :: targets(4)
type(normal_distribution) :: prior
integer :: chain_index

prior%mean_value = 0.0d0
prior%standard_deviation = 1.0d0
do chain_index = 1, size(targets)
	! Add each parameter's prior independently for this chain.
	call targets(chain_index)%add_parameter(prior)
	call targets(chain_index)%add_parameter(prior)
end do
call sample_chains(samplers, xdata, ydata, proposals, targets, niter=10000)
! Retrieve each result with samplers(chain_index)%get_chain().
```

The parameter count must match the custom model. Each proposal adapts its own
scale. As with `sample`, repeated calls append samples and initialize a new run
from the prior means; they do not resume from the last stored state. Omitting
`niter` uses 10,000 iterations per chain. Arrays of polymorphic objects must
have a common dynamic type within each array.

Likelihood evaluation uses a fused SIMD residual reduction without temporary
residual or log-density arrays. Outside an existing OpenMP parallel region, it
uses worker threads when the observation count reaches the target's public
`likelihood_parallel_threshold` (default 10,000). This default is a starting
point, not a benchmark-derived optimum: tune it for the workload, set it to zero
to request threading at every size, or to `huge(1_int32)` to disable threading
for ordinary dataset sizes. Nested worker teams are suppressed so independent
chains do not also launch likelihood teams. Reduction order can change floating
point results and, consequently, sampled trajectories.

Set `OMP_NUM_THREADS` to control concurrency, and avoid oversubscription from
threaded BLAS or custom models. Parallel chains require a compiler runtime with
thread-safe intrinsic `random_number`, such as GNU Fortran. Random-stream
assignment and reproducibility depend on that runtime and thread scheduling;
this API does not provide per-chain seeds or thread-count-independent results.
Without OpenMP enabled, the same API runs serially.

## Documentation
The generated API documentation is available [here](https://jchristopherson.github.io/fstats/).

## External Libraries
FSTATS uses [LINALG](https://github.com/jchristopherson/linalg) for linear algebra and [COLLECTIONS](https://github.com/jchristopherson/collections) for collection types. An optimized BLAS and LAPACK installation is recommended for best performance when using LINALG.
