
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
- [`Base.indexin`](api.md#Base.indexin)
- [`Base.indexin`](interfaces.md#Base.indexin-Tuple{GenomicVectors.AbstractGenomicVector,GenomicVectors.AbstractGenomicVector})
- [`Base.setdiff`](interfaces.md#Base.setdiff-Tuple{GenomicVectors.AbstractGenomicVector,GenomicVectors.AbstractGenomicVector})
- [`GenomicVectors.chromosomes`](api.md#GenomicVectors.chromosomes)
- [`GenomicVectors.chrpos`](api.md#GenomicVectors.chrpos)
- [`GenomicVectors.findoverlaps`](interfaces.md#GenomicVectors.findoverlaps)
- [`GenomicVectors.genoends`](interfaces.md#GenomicVectors.genoends)
- [`GenomicVectors.genopos`](api.md#GenomicVectors.genopos)
- [`GenomicVectors.genostarts`](interfaces.md#GenomicVectors.genostarts)
- [`GenomicVectors.nearest`](api.md#GenomicVectors.nearest)
- [`GenomicVectors.strands`](interfaces.md#GenomicVectors.strands)
- [`RLEVectors.ends`](api.md#RLEVectors.ends)
- [`RLEVectors.starts`](api.md#RLEVectors.starts)
- [`RLEVectors.widths`](api.md#RLEVectors.widths)


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


<a target='_blank' href='https://github.com/phaverty/GenomicVectors.jl/tree/5854804652d35a0ce3a952da821d6bab9bbcf80e/src/GenomeInfo_type.jl#L5-L23' class='documenter-source'>source</a><br>

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


<a target='_blank' href='https://github.com/phaverty/GenomicVectors.jl/tree/5854804652d35a0ce3a952da821d6bab9bbcf80e/src/GenomicPositions_type.jl#L5-L27' class='documenter-source'>source</a><br>

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


<a target='_blank' href='https://github.com/phaverty/GenomicVectors.jl/tree/5854804652d35a0ce3a952da821d6bab9bbcf80e/src/GenomicRanges_type.jl#L5-L35' class='documenter-source'>source</a><br>


<a id='Interfaces-1'></a>

## Interfaces

<a id='GenomicVectors.AbstractGenomicVector' href='#GenomicVectors.AbstractGenomicVector'>#</a>
**`GenomicVectors.AbstractGenomicVector`** &mdash; *Type*.



An AbstractGenomicVector is a Vector that describes positions or ranges in a single genome, in an arbitrary order. An AbstractGenomicVector must implement the GenomeInfo and GenoPos Interfaces. Sorting is by chromosome then by nucleotide position.


<a target='_blank' href='https://github.com/phaverty/GenomicVectors.jl/tree/5854804652d35a0ce3a952da821d6bab9bbcf80e/src/AbstractGenomicVector_type.jl#L7-L12' class='documenter-source'>source</a><br>


<a id='Accessing-position-info-1'></a>

## Accessing position info

<a id='RLEVectors.starts' href='#RLEVectors.starts'>#</a>
**`RLEVectors.starts`** &mdash; *Function*.



**RLEVectors**

`RLEVectors` is an alternate implementation of the Rle type from Bioconductor's IRanges package by H. Pages, P. Aboyoun and M. Lawrence. RLEVectors represent a vector with repeated values as the ordered set of values and repeat extents. In the field of genomics, data of various types measured across the ~3 billion letters in the human genome can often be represented in a few thousand runs. It is useful to know the bounds of genome regions covered by these runs, the values associated with these runs, and to be able to perform various mathematical operations on these values.

`RLEVectors` can be created from a single vector or a vector of values and a vector of run ends. In either case runs of values or zero length runs will be compressed out. RLEVectors can be expanded to a full vector with `collect`.

**Aliases**

Several aliases are defined for specific types of RLEVector (or collections thereof).

```
FloatRle              RLEVector{Float64,UInt32}
IntegerRle            RLEVector{Int64,UInt32}
BoolRle               RLEVector{Bool,UInt32}
StringRle             RLEVector{String,UInt32}
RLEVectorList{T1,T2}  Vector{ RLEVector{T1,T2} }
```

**Constructors**

`RLEVector`s can be created by specifying a vector to compress or the runvalues and run ends.

```
x = RLEVector([1,1,2,2,3,3,4,4,4])
x = RLEVector([4,5,6],[3,6,9])
```

**Describing `RLEVector` objects**

`RLEVector`s implement the usual descriptive functions for an array as well as some that are specific to the type.

  * `length(x)` The full length of the vector, uncompressed
  * `size(x)` Same as `length`, as for any other vector
  * `size(x,dim)` Returns `(length(x),1) for dim == 1`
  * `starts(x)` The index of the beginning of each run
  * `widths(x)` The width of each run
  * `ends(x)` The index of the end of each run
  * `values(x)` The data value for each run
  * `isempty(x)` Returns boolean, as for any other vector
  * `nrun(x)` Returns the number of runs represented in the array
  * `eltype(x)` Returns the element type of the runs
  * `endtype(x)` Returns the element type of the run ends


<a target='_blank' href='https://github.com/phaverty/RLEVectors.jl/tree/b4cd20186498d3e7709fb0044888cc6adad97388/src/describe.jl#L79' class='documenter-source'>source</a><br>


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
each(x)
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


<a target='_blank' href='https://github.com/phaverty/GenomicVectors.jl/tree/5854804652d35a0ce3a952da821d6bab9bbcf80e/src/utils.jl#L5-L51' class='documenter-source'>source</a><br>

<a id='RLEVectors.widths' href='#RLEVectors.widths'>#</a>
**`RLEVectors.widths`** &mdash; *Function*.



**RLEVectors**

`RLEVectors` is an alternate implementation of the Rle type from Bioconductor's IRanges package by H. Pages, P. Aboyoun and M. Lawrence. RLEVectors represent a vector with repeated values as the ordered set of values and repeat extents. In the field of genomics, data of various types measured across the ~3 billion letters in the human genome can often be represented in a few thousand runs. It is useful to know the bounds of genome regions covered by these runs, the values associated with these runs, and to be able to perform various mathematical operations on these values.

`RLEVectors` can be created from a single vector or a vector of values and a vector of run ends. In either case runs of values or zero length runs will be compressed out. RLEVectors can be expanded to a full vector with `collect`.

**Aliases**

Several aliases are defined for specific types of RLEVector (or collections thereof).

```
FloatRle              RLEVector{Float64,UInt32}
IntegerRle            RLEVector{Int64,UInt32}
BoolRle               RLEVector{Bool,UInt32}
StringRle             RLEVector{String,UInt32}
RLEVectorList{T1,T2}  Vector{ RLEVector{T1,T2} }
```

**Constructors**

`RLEVector`s can be created by specifying a vector to compress or the runvalues and run ends.

```
x = RLEVector([1,1,2,2,3,3,4,4,4])
x = RLEVector([4,5,6],[3,6,9])
```

**Describing `RLEVector` objects**

`RLEVector`s implement the usual descriptive functions for an array as well as some that are specific to the type.

  * `length(x)` The full length of the vector, uncompressed
  * `size(x)` Same as `length`, as for any other vector
  * `size(x,dim)` Returns `(length(x),1) for dim == 1`
  * `starts(x)` The index of the beginning of each run
  * `widths(x)` The width of each run
  * `ends(x)` The index of the end of each run
  * `values(x)` The data value for each run
  * `isempty(x)` Returns boolean, as for any other vector
  * `nrun(x)` Returns the number of runs represented in the array
  * `eltype(x)` Returns the element type of the runs
  * `endtype(x)` Returns the element type of the run ends


<a target='_blank' href='https://github.com/phaverty/RLEVectors.jl/tree/b4cd20186498d3e7709fb0044888cc6adad97388/src/describe.jl#L79' class='documenter-source'>source</a><br>


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
each(x)
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


<a target='_blank' href='https://github.com/phaverty/GenomicVectors.jl/tree/5854804652d35a0ce3a952da821d6bab9bbcf80e/src/utils.jl#L5-L51' class='documenter-source'>source</a><br>

<a id='RLEVectors.ends' href='#RLEVectors.ends'>#</a>
**`RLEVectors.ends`** &mdash; *Function*.



**RLEVectors**

`RLEVectors` is an alternate implementation of the Rle type from Bioconductor's IRanges package by H. Pages, P. Aboyoun and M. Lawrence. RLEVectors represent a vector with repeated values as the ordered set of values and repeat extents. In the field of genomics, data of various types measured across the ~3 billion letters in the human genome can often be represented in a few thousand runs. It is useful to know the bounds of genome regions covered by these runs, the values associated with these runs, and to be able to perform various mathematical operations on these values.

`RLEVectors` can be created from a single vector or a vector of values and a vector of run ends. In either case runs of values or zero length runs will be compressed out. RLEVectors can be expanded to a full vector with `collect`.

**Aliases**

Several aliases are defined for specific types of RLEVector (or collections thereof).

```
FloatRle              RLEVector{Float64,UInt32}
IntegerRle            RLEVector{Int64,UInt32}
BoolRle               RLEVector{Bool,UInt32}
StringRle             RLEVector{String,UInt32}
RLEVectorList{T1,T2}  Vector{ RLEVector{T1,T2} }
```

**Constructors**

`RLEVector`s can be created by specifying a vector to compress or the runvalues and run ends.

```
x = RLEVector([1,1,2,2,3,3,4,4,4])
x = RLEVector([4,5,6],[3,6,9])
```

**Describing `RLEVector` objects**

`RLEVector`s implement the usual descriptive functions for an array as well as some that are specific to the type.

  * `length(x)` The full length of the vector, uncompressed
  * `size(x)` Same as `length`, as for any other vector
  * `size(x,dim)` Returns `(length(x),1) for dim == 1`
  * `starts(x)` The index of the beginning of each run
  * `widths(x)` The width of each run
  * `ends(x)` The index of the end of each run
  * `values(x)` The data value for each run
  * `isempty(x)` Returns boolean, as for any other vector
  * `nrun(x)` Returns the number of runs represented in the array
  * `eltype(x)` Returns the element type of the runs
  * `endtype(x)` Returns the element type of the run ends


<a target='_blank' href='https://github.com/phaverty/RLEVectors.jl/tree/b4cd20186498d3e7709fb0044888cc6adad97388/src/describe.jl#L79' class='documenter-source'>source</a><br>


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
each(x)
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


<a target='_blank' href='https://github.com/phaverty/GenomicVectors.jl/tree/5854804652d35a0ce3a952da821d6bab9bbcf80e/src/utils.jl#L5-L51' class='documenter-source'>source</a><br>

<a id='GenomicVectors.genopos' href='#GenomicVectors.genopos'>#</a>
**`GenomicVectors.genopos`** &mdash; *Function*.



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
each(x)
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


<a target='_blank' href='https://github.com/phaverty/GenomicVectors.jl/tree/5854804652d35a0ce3a952da821d6bab9bbcf80e/src/utils.jl#L5-L51' class='documenter-source'>source</a><br>

<a id='GenomicVectors.chrpos' href='#GenomicVectors.chrpos'>#</a>
**`GenomicVectors.chrpos`** &mdash; *Function*.



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
each(x)
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


<a target='_blank' href='https://github.com/phaverty/GenomicVectors.jl/tree/5854804652d35a0ce3a952da821d6bab9bbcf80e/src/utils.jl#L5-L51' class='documenter-source'>source</a><br>

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
each(x)
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


<a target='_blank' href='https://github.com/phaverty/GenomicVectors.jl/tree/5854804652d35a0ce3a952da821d6bab9bbcf80e/src/utils.jl#L5-L51' class='documenter-source'>source</a><br>


<a id='Modifying-positions-1'></a>

## Modifying positions


```
slide
slide!
```


<a id='Querying-positions-1'></a>

## Querying positions


As in Bioconductor, location query operations discriminate between exact and overlapping matches. In addition to exact versus overlapping coordinates, exact matching includes strand, while overlap matching does not. In `GenomcVectors.jl`, the standard set operations use exact matching and custom overlap functions are defined for `AbstractGenomicVector`.


```
indexin
findin
overlap
overlaps # Hmm
hasoverlap
overlapin
overlapindex
nearest
```

