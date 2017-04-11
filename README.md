# GenomicVectors

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://phaverty.github.io/GenomicVectors.jl/stable)
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://phaverty.github.io/GenomicVectors.jl/latest)

`GenomicVectors` is an alternate implementation of the `GPos` and `GenomicRanges` types from
Bioconductor's GenomicRanges package by P. Aboyoun, H. Pages and M. Lawrence. These `GenomicPositions` and
`GenomicRanges` types are `Vectors` that serve as markers of locations on a given genome. They can be used
independently or as indices and/or annotation on other objects, such as the `GenomicDataFrame` implemented
here.

