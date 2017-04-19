# The GenomicVectors Types and Methods

## Index

```@index
```

## Types
```@docs
GenomeInfo
GenomicPositions
GenomicRanges
```

## Genome Location API
`GenomicVectors.jl` has ...
All `AbstractGenomicVector`s implement the API for `GenomeInfo` for access to their genome descriptions.


## Accessing position info
```@docs
GenomicVectors.starts
GenomicVectors.widths
GenomicVectors.ends
genopos
chrpos
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
does not. In `GenomcVectors.jl`, the standard set operations use exact matching and custom overlap
functions are defined for `AbstractGenomicVector`.

### Overlap function
```@docs
overlap
overlapin
overlapindex
nearest
```
