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

## A round trip to R
```julia
using RCall
chrinfo = GenomeInfo("hg19",["chr1","chr2","chrX"],Int64[3e5,2e5,1e4])
chrs = ["chr1","chr2","chr2","chrX"]
s = [100, 200, 300, 400]
e = [120, 240, 350, 455]
x = GenomicRanges(chrs,s,e,chrinfo)
y = RObject(x)
@rput y
R"z = y + 2L"
@rget z
```
