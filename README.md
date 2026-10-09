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

## Documentation
The generated API documentation is available [here](https://jchristopherson.github.io/fstats/).

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

Likelihood evaluation uses scaled, two-pass SIMD residual reductions without temporary
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

## Distribution Modules
Each concrete distribution has its own module and source file in
`src/distributions`. The module name is `fstats_` followed by the type name,
for example:

```fortran
use fstats_distributions, only : distribution
use fstats_normal_distribution, only : normal_distribution
```

`fstats_distributions` contains the abstract `distribution` and
`multivariate_distribution` types, their interfaces, and shared supporting
routines. It no longer exports concrete distribution types. Replace direct
imports of concrete types from that module with their individual modules,
or continue to use `use fstats` to access all distributions unchanged.

Each built-in distribution implements its own log-density and tail bindings;
the abstract bases supply generic compatibility fallbacks for custom laws.

## Numerical Behavior
Distributions expose `log_pdf`, `survival`, `log_cdf`, and `log_survival` in
addition to `pdf` and `cdf`. Prefer log densities for posterior calculations
and direct survival functions for small upper-tail probabilities. Built-in
laws use analytic log densities and direct beta/gamma or normal tails. Discrete
PDFs are probability masses at integer-valued points; their CDFs include all
masses at or below the supplied point, and survival means `P(X > x)`.

Density and tail evaluations with invalid parameters or NaN inputs return NaN. Outside
support, PDFs return zero and log densities return negative infinity. A
support-endpoint density singularity can legitimately return positive infinity.
Custom distributions have compatibility fallbacks that take logs of their PDFs
or CDFs; override these methods for accurate extreme tails. Subclasses that
change a built-in probability law must also override its log-density/tail methods.

`log_regularized_beta`, `log_regularized_gamma_lower`, and
`log_regularized_gamma_upper` avoid underflow in special-function probabilities.
The regularized routines accept an optional `max_iterations` limit and return
NaN on invalid input or failure to converge. `stable_log1p` and `stable_expm1`
retain small differences lost by `log(1 + x)` and `exp(x) - 1`.

MCMC stores log variance and evaluates the default Gaussian likelihood without
exponentiating it. `likelihood` still accepts variance, but nonpositive or
nonfinite variance now returns NaN rather than being silently clamped.
`likelihood_log_variance` accepts log variance directly. `log_posterior` uses
`evaluate_log_variance_prior`, whose density includes the variance-transform
Jacobian. Custom targets changing the likelihood or variance prior must override
these log-coordinate bindings; overriding only the variance-space methods no
longer changes sampling. Initial states must have finite log posterior; invalid
or impossible proposals are rejected. Correcting the Jacobian changes sampled
variance distributions relative to earlier versions.

Descriptive statistics use scaled/compensated reductions. Empty statistics,
sample moments with fewer than two observations, and nonfinite observations
return NaN; constant-data correlation also returns NaN. Quantiles require
`0 <= q <= 1`, trimming requires `0 <= p < 0.5`, and pooled variance requires
at least two observations per group. These contracts replace some earlier zero
results or silent clamping. Standard deviation can remain finite even when
variance overflows; truly unrepresentable moments return infinity.

Multivariate-normal log densities use Cholesky solves and log determinants.
Regression covariance uses an SVD of the design matrix, not its normal equations.
`regression_covariance` optionally reports numerical rank and condition number;
nonlinear fits expose the same diagnostics in `convergence_info` when covariance
is requested. Rank-deficient covariance is a pseudocovariance: discarded
directions are not evidence of zero parameter uncertainty. Nonlinear covariance
requires additional evaluations of the final Jacobian, counted in the reported
function evaluations. Ordinary regression coefficients retain their existing
parameterization; scaling polynomial predictors remains advisable.

## Experimental Design Modules
All experimental design types are declared in `fstats_experimental_design`.
Procedures are organized into the following modules in `src/design`:

- `fstats_doe_designs`: design generation, size calculations, and efficiency assessment.
- `fstats_doe_coding`: natural-to-coded and coded-to-natural factor conversions.
- `fstats_doe_models`: model fitting and evaluation.
- `fstats_doe_diagnostics`: goodness of fit, residual analysis, model comparison, and model ANOVA.
- `fstats_doe_prediction`: predictions and uncertainty intervals.
- `fstats_doe_response_surface`: response-surface optimization.

`use fstats` continues to expose all these types and procedures. Code that
previously imported procedures directly from `fstats_experimental_design` must
instead import them from the corresponding module above, or use `fstats`.
Import shared types from `fstats_experimental_design` when using individual
procedure modules. This organization does not change numerical behavior;
existing approximations and placeholder calculations remain unchanged.

## External Libraries
FSTATS uses [LINALG](https://github.com/jchristopherson/linalg) for linear algebra and [COLLECTIONS](https://github.com/jchristopherson/collections) for collection types. An optimized BLAS and LAPACK installation is recommended for best performance when using LINALG.
