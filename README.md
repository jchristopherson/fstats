# fstats
FSTATS is a modern Fortran statistical library containing routines for computing basic statistical properties, hypothesis testing, regression, special functions, and even experimental design.

## Status
[![CMake](https://github.com/jchristopherson/fstats/actions/workflows/cmake.yml/badge.svg)](https://github.com/jchristopherson/fstats/actions/workflows/cmake.yml)
[![Actions Status](https://github.com/jchristopherson/fstats/workflows/fpm/badge.svg)](https://github.com/jchristopherson/fstats/actions)

## Building FSTATS
[CMake](https://cmake.org/)This library can be built using CMake.  For instructions see [Running CMake](https://cmake.org/runningcmake/).

[FPM](https://github.com/fortran-lang/fpm) can also be used to build this library using the provided fpm.toml.
```txt
fpm build
```
The FSTATS library can be used within your FPM project by adding the following to your fpm.toml file.
```toml
[dependencies]
fstats = { git = "https://github.com/jchristopherson/fstats" }
```

## Documentation
Documentation can be found [here](https://jchristopherson.github.io/fstats/).

## External Libraries
Here is a list of external code libraries utilized by this library.  The CMake build script will include these dependencies automatically; however, it is highly recommended that an optimized BLAS and LAPACK already reside on your system for best performance (used by LINALG for linear algebra calculations).
- [FERROR](https://github.com/jchristopherson/ferror)
- [LINALG](https://github.com/jchristopherson/linalg)
- [COLLECTIONS](https://github.com/jchristopherson/collections)

## Capabilities
- Descriptive Statistics
- Hypothesis Testing
- ANOVA
- Regression
- Bootstrapping
- Sampling
- Special Functions
- Allan Variance
- Smoothing
- Experimental Design
- Markov Chain Monte Carlo
- Interpolation
