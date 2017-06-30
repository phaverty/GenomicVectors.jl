# The GenomicVectors API

```@contents
```

## Types
```@docs
AbstractGenomicVector
GenomeInfo
GenomicPositions
GenomicRanges
```

## Interfaces and Accessing position info
```@docs
chromosomes
```

## Modifying positions
```@docs
slide
slide!
```

## Querying positions
As in Bioconductor, location query operations discriminate between exact and overlapping matches. In
addition to exact versus overlapping coordinates, exact matching includes strand, while overlap matching
does not. In `GenomicVectors.jl`, the standard set operations use exact matching and custom overlap functions are defined for `AbstractGenomicVector`.

```@docs
findoverlaps
```
