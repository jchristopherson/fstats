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

## Documentation
The generated API documentation is available [here](https://jchristopherson.github.io/fstats/).

## External Libraries
FSTATS uses [LINALG](https://github.com/jchristopherson/linalg) for linear algebra and [COLLECTIONS](https://github.com/jchristopherson/collections) for collection types. An optimized BLAS and LAPACK installation is recommended for best performance when using LINALG.
