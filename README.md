# GenomicVectors

[![Build Status](https://travis-ci.org/phaverty/GenomicVectors.jl.svg?branch=master)](https://travis-ci.org/phaverty/GenomicVectors.jl)
[![Coverage Status](https://codecov.io/github/phaverty/GenomicVectors.jl/coverage.svg?branch=master)](https://codecov.io/github/phaverty/GenomicVectors.jl?branch=master)
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://phaverty.github.io/RLEVectors.jl/latest)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://phaverty.github.io/RLEVectors.jl/stable)

`GenomicVectors` is an alternate implementation of the `GPos` and `GenomicRanges` types from
Bioconductor's GenomicRanges package by P. Aboyoun, H. Pages and M. Lawrence. These `GenomicPositions` and
`GenomicRanges` types are `Vectors` that serve as markers of locations on a given genome. They can be used
independently or as indices and/or annotation on other objects.

## Installation

Releases of GenomicVectors version 1.0 and above are registered and made available
to install through BioJulia's package registry. Julia's package manager only
monitors the "General" package repository by default. So before you start, you
should tell julia about the existence of the BioJulia package registry.

Start a julia terminal, hit the `]` key to enter pkg mode (you should see the
prompt change from `julia>` to `pkg>` ), then enter the following command:

```julia
registry add https://github.com/BioJulia/BioJuliaRegistry.git
```

After you've added the registry, you can install BioSequences from the julia
REPL. Press `]` to enter pkg mode again, and enter the following:

```julia
add GenomicVectors
```

If you are interested in the cutting edge of the development, please check out
the master branch to try new features before release.
