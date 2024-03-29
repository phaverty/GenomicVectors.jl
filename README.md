# GenomicVectors

[![CI Testing](https://github.com/phaverty/GenomicVectors.jl/workflows/CI/badge.svg)](https://github.com/phaverty/GenomicVectors.jl/actions?query=workflow%3ACI+branch%3Amain)
[![Coverage Status](https://codecov.io/github/phaverty/GenomicVectors.jl/coverage.svg?branch=master)](https://codecov.io/github/phaverty/GenomicVectors.jl?branch=master)
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://phaverty.github.io/GenomicVectors.jl/latest)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://phaverty.github.io/GenomicVectors.jl/stable)

`GenomicVectors` is an alternate implementation of the `GPos` and `GenomicRanges` types from
Bioconductor's GenomicRanges package by P. Aboyoun, H. Pages and M. Lawrence. These `GenomicPositions` and
`GenomicRanges` types are `Vectors` that serve as markers of locations on a given genome. They can be used
independently or as indices and/or annotation on other objects.
