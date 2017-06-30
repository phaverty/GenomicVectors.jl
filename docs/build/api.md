
<a id='The-GenomicVectors-API-1'></a>

# The GenomicVectors API

- [The GenomicVectors API](api.md#The-GenomicVectors-API-1)
    - [Types](api.md#Types-1)
    - [Interfaces and Accessing position info](api.md#Interfaces-and-Accessing-position-info-1)
    - [Modifying positions](api.md#Modifying-positions-1)
    - [Querying positions](api.md#Querying-positions-1)
    - [Interfaces](interfaces.md#Interfaces-1)
- [Ideas for future development](future.md#Ideas-for-future-development-1)
    - [Maintaining performance via sortedness and views](future.md#Maintaining-performance-via-sortedness-and-views-1)
- [TODO](TODO.md#TODO-1)
    - [Features](TODO.md#Features-1)
    - [Decisions](TODO.md#Decisions-1)
    - [Improvements](TODO.md#Improvements-1)
    - [Bugs](TODO.md#Bugs-1)
- [Release Notes](NEWS.md#Release-Notes-1)
    - [Version 0.0.1: Initial Public Release](NEWS.md#Version-0.0.1:-Initial-Public-Release-1)
    - [Version 0.0.2: changes to indexing a GenomicRanges](NEWS.md#Version-0.0.2:-changes-to-indexing-a-GenomicRanges-1)
    - [Version 0.0.3: Update to scalar indexing, behaves like Vector{Interval}](NEWS.md#Version-0.0.3:-Update-to-scalar-indexing,-behaves-like-Vector{Interval}-1)
    - [Version 0.1.0: Introducing `GenomicTable`](NEWS.md#Version-0.1.0:-Introducing-GenomicTable-1)
- [GenomicVectors](index.md#GenomicVectors-1)
    - [Introduction](index.md#Introduction-1)
    - [Implementation](index.md#Implementation-1)
    - [Working with locations](index.md#Working-with-locations-1)
    - [Ordering](index.md#Ordering-1)
    - [Intersection / overlap operations](index.md#Intersection-/-overlap-operations-1)


<a id='Types-1'></a>

## Types

<a id='GenomicVectors.AbstractGenomicVector' href='#GenomicVectors.AbstractGenomicVector'>#</a>
**`GenomicVectors.AbstractGenomicVector`** &mdash; *Type*.



An AbstractGenomicVector is a Vector that describes positions or ranges in a single genome, in an arbitrary order. An AbstractGenomicVector must implement the GenomeInfo and GenoPos Interfaces. Sorting is by chromosome then by nucleotide position.


<a target='_blank' href='https://github.com/phaverty/GenomicVectors.jl/tree/bc4d311571f27c868b680e986018f31466704dcc/src/AbstractGenomicVector_type.jl#L5-L10' class='documenter-source'>source</a><br>

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


<a target='_blank' href='https://github.com/phaverty/GenomicVectors.jl/tree/bc4d311571f27c868b680e986018f31466704dcc/src/GenomeInfo_type.jl#L5-L23' class='documenter-source'>source</a><br>

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


<a target='_blank' href='https://github.com/phaverty/GenomicVectors.jl/tree/bc4d311571f27c868b680e986018f31466704dcc/src/GenomicPositions_type.jl#L5-L27' class='documenter-source'>source</a><br>

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

Getting/setting by a scalar gives/takes a GenomicFeatures.Interval. The leftposition and rightposition in this Interval must be in genome location units and correspond to the same chromosome. The seqname must match the genome of the GenomicRanges. Outgoing Intervals will have the index `i` as their metadata. This makes it possible to obtain the original ordering if Intervals after conversion to, say, an IntervalCollection. Any metadata  for an incoming Interval is ignored.

The `each` function produces an iterator of (start,end) two-tuples in genome location units. This is use for many internal functions, like sorting. This is intentionally similar to `RLEVectors.each`.


<a target='_blank' href='https://github.com/phaverty/GenomicVectors.jl/tree/bc4d311571f27c868b680e986018f31466704dcc/src/GenomicRanges_type.jl#L5-L35' class='documenter-source'>source</a><br>


<a id='Interfaces-and-Accessing-position-info-1'></a>

## Interfaces and Accessing position info

<a id='GenomicVectors.chromosomes' href='#GenomicVectors.chromosomes'>#</a>
**`GenomicVectors.chromosomes`** &mdash; *Function*.



**The GenoPos Interface**

Provides access to positional information in the linearized genome or in chromosome coordinate (e.g. chr4:1000-1020). This interface requires a type to implement the a method on the non-copying accessors `_genostarts`, `_genoends` and `_strands` as well as the GenomeInfo Interface.

```
starts(x)
```

Get the starting nucleotide index for each range/position relative to the chromosome on which they lie.

```
ends(x)
```

Get the ending nucleotide index for each range/position relative to the chromosome on which they lie.

```
widths(x)
```

Get the distance, between the start and end nucleotide of the range, An `RleVector` of 1s for single-nucleotide positions.

```
chromosomes(x)
chromosomes(genopos, chrinfo)
```

Get the name of the chromosome for each range/position.

```
genostarts(x)
```

Get the starting nucleotide index for each range/position in the linearized genome.

```
genoends(x)
```

Get the ending nucleotide index for each range/position in the linearized genome.

```
strands(x)
```

Get the DNA strand for each range/position, pass by copy.

```
eachrange(x)
```

Return an iterator that returns tuples of genostart and genoend pairs.

**Utility Functions**

Much of this functionality is derived from a few utility functions:

```
chrpos(genopos, chrinfo)
```

Given positions in the linear genome, calculate the position on the relevant chromosome.

```
chrindex(genopos, chrinfo)
```

Given positions in the linear genome, determine the corresponding chromosomes and return the indices of the chromosome in chrinfo.

```
genopos(chrpos, chromosomes, chrinfo)
```

Given chromosome and chromosome position information and a description of the chromosomes (a GenoPos object), calculate the corresponding positions in the linear genome.


<a target='_blank' href='https://github.com/phaverty/GenomicVectors.jl/tree/bc4d311571f27c868b680e986018f31466704dcc/src/utils.jl#L5-L51' class='documenter-source'>source</a><br>


<a id='Modifying-positions-1'></a>

## Modifying positions


```
slide
slide!
```


<a id='Querying-positions-1'></a>

## Querying positions


As in Bioconductor, location query operations discriminate between exact and overlapping matches. In addition to exact versus overlapping coordinates, exact matching includes strand, while overlap matching does not. In `GenomicVectors.jl`, the standard set operations use exact matching and custom overlap functions are defined for `AbstractGenomicVector`.

<a id='GenomicVectors.findoverlaps' href='#GenomicVectors.findoverlaps'>#</a>
**`GenomicVectors.findoverlaps`** &mdash; *Function*.



**Searching GenomicVectors**

Matching functions in `GenomicVectors` can perform overlap matching, rather than exact matching when given the extra argument `exact=false`. In either case, the genome strand is never considered. Other types of [overlaps](https://en.wikipedia.org/wiki/Allen%27s_interval_algebra) may be supported in the future.

```
findoverlaps(x::AbstractGenomicVector,y::AbstractGenomicVector)
```

Creates a `GenomicFeatures.IntersectIterator` from two `AbstractGenomicVectors`, much like the BioConductor `findOverlaps`. This is the kernel of the other search/set operations.

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


<a target='_blank' href='https://github.com/phaverty/GenomicVectors.jl/tree/bc4d311571f27c868b680e986018f31466704dcc/src/AbstractGenomicVector_type.jl#L68-L93' class='documenter-source'>source</a><br>

