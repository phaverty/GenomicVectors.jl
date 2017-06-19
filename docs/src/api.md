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

## Interfaces
```@docs
AbstractGenomicVector
```

## Accessing position info
```@docs
starts(AbstractGenomicVector)
widths(AbstractGenomicVector)
ends(AbstractGenomicVector)
chromosomes(AbstractGenomicVector)
genopos
chrpos
chrindex
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

```@docs
findoverlaps
indexin
findin
nearest
```
