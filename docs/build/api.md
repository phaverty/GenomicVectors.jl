
<a id='The-GenomicVectors-Types-and-Methods-1'></a>

# The GenomicVectors Types and Methods


<a id='Index-1'></a>

## Index

- [`GenomicVectors.AbstractGenomicVector`](api.md#GenomicVectors.AbstractGenomicVector)
- [`GenomicVectors.GenomeInfo`](api.md#GenomicVectors.GenomeInfo)
- [`GenomicVectors.GenomicPositions`](api.md#GenomicVectors.GenomicPositions)
- [`GenomicVectors.GenomicRanges`](api.md#GenomicVectors.GenomicRanges)
- [`Base.findin`](api.md#Base.findin)
- [`Base.findin`](interfaces.md#Base.findin-Tuple{GenomicVectors.AbstractGenomicVector,GenomicVectors.AbstractGenomicVector})
- [`Base.in`](interfaces.md#Base.in-Tuple{GenomicVectors.AbstractGenomicVector,GenomicVectors.AbstractGenomicVector})
- [`Base.indexin`](interfaces.md#Base.indexin-Tuple{GenomicVectors.AbstractGenomicVector,GenomicVectors.AbstractGenomicVector})
- [`Base.indexin`](api.md#Base.indexin)
- [`Base.setdiff`](interfaces.md#Base.setdiff-Tuple{GenomicVectors.AbstractGenomicVector,GenomicVectors.AbstractGenomicVector})
- [`GenomicVectors.chrindex`](api.md#GenomicVectors.chrindex)
- [`GenomicVectors.chromosomes`](interfaces.md#GenomicVectors.chromosomes)
- [`GenomicVectors.chrpos`](api.md#GenomicVectors.chrpos)
- [`GenomicVectors.findoverlaps`](api.md#GenomicVectors.findoverlaps)
- [`GenomicVectors.genoends`](interfaces.md#GenomicVectors.genoends)
- [`GenomicVectors.genopos`](api.md#GenomicVectors.genopos)
- [`GenomicVectors.genostarts`](interfaces.md#GenomicVectors.genostarts)
- [`GenomicVectors.nearest`](api.md#GenomicVectors.nearest)
- [`GenomicVectors.strands`](interfaces.md#GenomicVectors.strands)


<a id='Types-1'></a>

## Types

<a id='GenomicVectors.GenomeInfo' href='#GenomicVectors.GenomeInfo'>#</a>
**`GenomicVectors.GenomeInfo`** &mdash; *Type*.



**GenomeInfo Type**

Describes a genome as a collection of chromosomes arranged end-to-end in a specific order such that the index of a nucleotide in any chromosome may be described by as single integer.

Indexing a `GenomeInfo` returns the genopos end of the indexed chromosome.

**Examples**

```julia
chrinfo = GenomeInfo("hg19",["chr1","chr2","chrX"],Int64[3e5,2e5,1e4])
genome(chrinfo)
chr_names(chrinfo)
chr_lengths(chrinfo)
chr_ends(chrinfo)
chr_offsets(chrinfo)
chrinfo[2] # 5e5
```


<a target='_blank' href='https://github.com/phaverty/GenomicVectors.jl/tree/897db59a49f6b4c76c8db11521f107638c4aac8a/src/GenomeInfo_type.jl#L5-L23' class='documenter-source'>source</a><br>

<a id='GenomicVectors.GenomicPositions' href='#GenomicVectors.GenomicPositions'>#</a>
**`GenomicVectors.GenomicPositions`** &mdash; *Type*.



**GenomicPositions Type**

Represents single-nucleotide positions in a genome.

This type uses its (immutable) `GenomeInfo` slot object to describe corresponding genome. Therefore, positions can be expressed relative to this concatenated, linearized genome or relative to the chromosome containing a given position.

By convention, all postions in a `GenomicPositions` are considered to be on the plus strand.

**Examples**

```julia
    genomeinfo = GenomeInfo("hg19",["chr1","chr2","chrX"],Int64[3e5,2e5,1e4])
    chrs = ["chr2","chr1","chr2","chrX"]
    pos = Int64[3e4,4.2e3,1.9e5,1e4]
    gpos = genopos(pos,chrs,chrinfo)
    x = GenomicPositions(pos,chrs,genomeinfo)
    y = GenomicPositions(gpos,genomeinfo)
    same_genome(x, y)
    convert(DataTable, y)
```


<a target='_blank' href='https://github.com/phaverty/GenomicVectors.jl/tree/897db59a49f6b4c76c8db11521f107638c4aac8a/src/GenomicPositions_type.jl#L5-L27' class='documenter-source'>source</a><br>

<a id='GenomicVectors.GenomicRanges' href='#GenomicVectors.GenomicRanges'>#</a>
**`GenomicVectors.GenomicRanges`** &mdash; *Type*.



**GenomicRanges Type**

Represents closed ranges in a genome. This type uses its (immutable) `GenomeInfo` slot object to describe corresponding genome and positions can be expressed relative to this concatenated, linearized genome or relative to the chromosome containing a given position.

**Examples**

```julia
    chrinfo = GenomeInfo("hg19",["chr1","chr2","chrX"],Int64[3e5,2e5,1e4])
    chrs = ["chr1","chr2","chr2","chrX"]
    starts = [100, 200, 300, 400]
    ends = [120, 240, 350, 455]
    gr = GenomicRanges(chrs,starts,ends,chrinfo)
```

**Indexing**

Indexing a `GenomicRanges` with an array produces a new `GenomicRanges`.

Getting/setting by a scalar gives/takes a Bio.Intervals.Interval. The leftposition and rightposition in this Interval must be in genome location units and correspond to the same chromosome. The seqname must match the genome of the GenomicRanges. Outgoing Intervals will have the index `i` as their metadata. This makes it possible to obtain the original ordering if Intervals after conversion to, say, an IntervalCollection. Any metadata  for an incoming Interval is ignored.

The `each` function produces an iterator of (start,end) two-tuples in genome location units. This is use for many internal functions, like sorting. This is intentionally similar to `RLEVectors.each`.


<a target='_blank' href='https://github.com/phaverty/GenomicVectors.jl/tree/897db59a49f6b4c76c8db11521f107638c4aac8a/src/GenomicRanges_type.jl#L5-L35' class='documenter-source'>source</a><br>


<a id='Interfaces-1'></a>

## Interfaces

<a id='GenomicVectors.AbstractGenomicVector' href='#GenomicVectors.AbstractGenomicVector'>#</a>
**`GenomicVectors.AbstractGenomicVector`** &mdash; *Type*.



An AbstractGenomicVector is a Vector that describes positions or ranges in a single genome, in an arbitrary order. An AbstractGenomicVector must implement the GenomeInfo and GenoPos Interfaces. Sorting is by chromosome then by nucleotide position.


<a target='_blank' href='https://github.com/phaverty/GenomicVectors.jl/tree/897db59a49f6b4c76c8db11521f107638c4aac8a/src/AbstractGenomicVector_type.jl#L7-L12' class='documenter-source'>source</a><br>


<a id='Accessing-position-info-1'></a>

## Accessing position info


```
starts(AbstractGenomicVector)
widths(AbstractGenomicVector)
ends(AbstractGenomicVector)
chromosomes(AbstractGenomicVector)
genopos
chrpos
chrindex
```


<a id='Modifying-positions-1'></a>

## Modifying positions


```
slide
slide!
```


<a id='Querying-positions-1'></a>

## Querying positions


As in Bioconductor, location query operations discriminate between exact and overlapping matches. In addition to exact versus overlapping coordinates, exact matching includes strand, while overlap matching does not. In `GenomcVectors.jl`, the standard set operations use exact matching and custom overlap functions are defined for `AbstractGenomicVector`.

<a id='GenomicVectors.findoverlaps' href='#GenomicVectors.findoverlaps'>#</a>
**`GenomicVectors.findoverlaps`** &mdash; *Function*.



**Searching GenomicVectors**

Matching functions in `GenomicVectors` can perform overlap matching, rather than exact matching when given the extra argument `exact=false`. In either case, the genome strand is never considered.

```
findoverlaps(x::AbstractGenomicVector,y::AbstractGenomicVector)
```

Creates a `Bio.Intervals.IntersectIterator` from two `AbstractGenomicVectors`, much like the BioConductor `findOverlaps`. `Bio.Intervals` calls this function `intersect`, but I would expect `intersect` to have the same behavior as base, returning a subset copy of the first argument. `findoverlaps` is the kernel of the other search/set operations.

```
findin(x::AbstractGenomicVector,y::AbstractGenomicVector,exact::Bool=true)

indexin(x::AbstractGenomicVector, y::AbstractGenomicVector, exact::Bool=true)

in(x::AbstractGenomicVector, y::AbstractGenomicVector, exact::Bool=true)
```

The `AbstractGenomicVector` method on `in` is vectorized and returns a `BitArray` that is `true` for each element of `x` that is in the set `y`.

```
intersect(x::AbstractGenomicVector, y::AbstractGenomicVector, exact::Bool=true)

setdiff(x::AbstractGenomicVector, y::AbstractGenomicVector, exact::Bool=true)
```


<a target='_blank' href='https://github.com/phaverty/GenomicVectors.jl/tree/897db59a49f6b4c76c8db11521f107638c4aac8a/src/AbstractGenomicVector_type.jl#L66-L92' class='documenter-source'>source</a><br>

<a id='Base.indexin' href='#Base.indexin'>#</a>
**`Base.indexin`** &mdash; *Function*.



**Searching GenomicVectors**

Matching functions in `GenomicVectors` can perform overlap matching, rather than exact matching when given the extra argument `exact=false`. In either case, the genome strand is never considered.

```
findoverlaps(x::AbstractGenomicVector,y::AbstractGenomicVector)
```

Creates a `Bio.Intervals.IntersectIterator` from two `AbstractGenomicVectors`, much like the BioConductor `findOverlaps`. `Bio.Intervals` calls this function `intersect`, but I would expect `intersect` to have the same behavior as base, returning a subset copy of the first argument. `findoverlaps` is the kernel of the other search/set operations.

```
findin(x::AbstractGenomicVector,y::AbstractGenomicVector,exact::Bool=true)

indexin(x::AbstractGenomicVector, y::AbstractGenomicVector, exact::Bool=true)

in(x::AbstractGenomicVector, y::AbstractGenomicVector, exact::Bool=true)
```

The `AbstractGenomicVector` method on `in` is vectorized and returns a `BitArray` that is `true` for each element of `x` that is in the set `y`.

```
intersect(x::AbstractGenomicVector, y::AbstractGenomicVector, exact::Bool=true)

setdiff(x::AbstractGenomicVector, y::AbstractGenomicVector, exact::Bool=true)
```


<a target='_blank' href='https://github.com/phaverty/GenomicVectors.jl/tree/897db59a49f6b4c76c8db11521f107638c4aac8a/src/AbstractGenomicVector_type.jl#L66-L92' class='documenter-source'>source</a><br>

<a id='Base.findin' href='#Base.findin'>#</a>
**`Base.findin`** &mdash; *Function*.



**Searching GenomicVectors**

Matching functions in `GenomicVectors` can perform overlap matching, rather than exact matching when given the extra argument `exact=false`. In either case, the genome strand is never considered.

```
findoverlaps(x::AbstractGenomicVector,y::AbstractGenomicVector)
```

Creates a `Bio.Intervals.IntersectIterator` from two `AbstractGenomicVectors`, much like the BioConductor `findOverlaps`. `Bio.Intervals` calls this function `intersect`, but I would expect `intersect` to have the same behavior as base, returning a subset copy of the first argument. `findoverlaps` is the kernel of the other search/set operations.

```
findin(x::AbstractGenomicVector,y::AbstractGenomicVector,exact::Bool=true)

indexin(x::AbstractGenomicVector, y::AbstractGenomicVector, exact::Bool=true)

in(x::AbstractGenomicVector, y::AbstractGenomicVector, exact::Bool=true)
```

The `AbstractGenomicVector` method on `in` is vectorized and returns a `BitArray` that is `true` for each element of `x` that is in the set `y`.

```
intersect(x::AbstractGenomicVector, y::AbstractGenomicVector, exact::Bool=true)

setdiff(x::AbstractGenomicVector, y::AbstractGenomicVector, exact::Bool=true)
```


<a target='_blank' href='https://github.com/phaverty/GenomicVectors.jl/tree/897db59a49f6b4c76c8db11521f107638c4aac8a/src/AbstractGenomicVector_type.jl#L66-L92' class='documenter-source'>source</a><br>

<a id='GenomicVectors.nearest' href='#GenomicVectors.nearest'>#</a>
**`GenomicVectors.nearest`** &mdash; *Function*.



```
function nearest(query::GenomicPositions, target::GenomicPositions)
```

For each `query` finds index in `target` that is nearest on the same chromosome. If no match on the same chromosome exists, the index will be 0.


<a target='_blank' href='https://github.com/phaverty/GenomicVectors.jl/tree/897db59a49f6b4c76c8db11521f107638c4aac8a/src/GenomicPositions_type.jl#L148-L152' class='documenter-source'>source</a><br>

