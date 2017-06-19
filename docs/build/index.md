
<a id='GenomicVectors-1'></a>

# GenomicVectors


<a id='Introduction-1'></a>

## Introduction


`GenomicVectors` is an alternate implementation of the `GPos` and `GenomicRanges` types from Bioconductor's GenomicRanges package by P. Aboyoun, H. Pages and M. Lawrence. These `GenomicPositions` and `GenomicRanges` types are `Vectors` that serve as markers of locations on a given genome. They can be used independently or as indices and/or annotation on other objects.


These types are commonly used in conjuction with `RLEVectors` from the package of the same name, which also often contain data arrayed along a genome. For example `GenomicDataFrame` may use a `GPos` or `GenomicRanges` as a row index and use `RLEVector` objects for data columns. This is a common method of storing segmented DNA copy number for multiple samples in R's `genoset` package.


<a id='Implementation-1'></a>

## Implementation


These `Vector` types each contain a `GenomicInfo` object, which annotates the names of the relevant genome and its chromosomes as well as the lengths of each chromosome. Operations on two or more of these genomic vectors require that they contain identical `GenomeInfo` objects.


The primary "innovation" of these types is that genome locations are stored as the 1-based index into the linear genome of concatenated chromosomes described by the immutable `GenomeInfo` object. The relevant chromosome for this `genopos` can be looked up efficiently as the `GenomeInfo` holds the `cumsum` of the chromosome lengths. Using binary search and the usual optimization for a sorted vector of queries, these lookups are O(numpositions * log(numchromosomes)) in the worst case and O(numpositions) in the best case.


<a id='Creation-1'></a>

### Creation


`GenomeInfo` and `GenomicPositions` objects can be created as follows:


```julia
using GenomicVectors

chrinfo = GenomeInfo("hg19",["chr1","chr2","chrX"],Int64[3e5,2e5,1e4])
chrs = ["chr1","chr2","chr2","chrX"]
pos = Int64[3e5,1.8e5,1.9e5,1e4]

x = GenomicPositions(pos,chrs,chrinfo)

gpos = genopos(pos,chrs,chrinfo)
y = GenomicPositions(gpos,chrinfo)
```


<a id='Describing-1'></a>

### Describing


`GenomeInfo` objects have various accessors that describe the chomosomes and their boundaries in the concatenated, linear genome. They are also immutable so that the meaning of indices into a concatenated genome cannot change unexpectedly. Objects that contain a `GenomeInfo` implement the same methods.


```julia
chrinfo = GenomeInfo("hg19",["chr1","chr2","chrX"],Int64[3e5,2e5,1e4])
genome(chrinfo)
chr_names(chrinfo)
chr_ends(chrinfo)
chr_lengths(chrinfo)
chr_offsets(chrinfo)

chrs = ["chr1","chr2","chr2","chrX"]
pos = Int64[3e5,1.8e5,1.9e5,1e4]
x = GenomicPositions(pos,chrs,chrinfo)
chr_names(x)
chr_ends(x)
chr_lengths(x)
chr_offsets(x)
```


Similarly, `GenomicPositions` and similar objects can describe the positions they represent. `GenomicRanges` and `GenomicPositions` have the same API, for interchangeability, even when functions are more naturally position- or range-related. Some of these functions are shared with `RLEVectors` as ranges and runs have much in common.


```julia
chrinfo = GenomeInfo("hg19",["chr1","chr2","chrX"],Int64[3e5,2e5,1e4])
chrs = ["chr1","chr2","chr2","chrX"]
pos = Int64[3e5,1.8e5,1.9e5,1e4]
x = GenomicPositions(pos,chrs,chrinfo)
chrpos(x)
starts(x)
ends(x)
widths(x)
```


<a id='Working-with-locations-1'></a>

## Working with locations


Genome locations may be modified (slid along the genome, resized, etc.) in various ways, but operations that would exceed the bounds of the relevant chromosome generate errors.


Queries of a genomic vector are only relevant in the context of a specific genome. Therefore, comparisons are only implemented between the collection of types that can be checked for identical `GenomeInfo` information.


Standard vector item matching operations require identity of items, rather than overlapping ranges. Overlap operations are implemented as distinct functions, as in Bioconductor.


```julia
chrinfo = GenomeInfo("GRCh38",["chr1","chr2","chr3","chrX","chrY","chrM"],Int64[3e5,2e5,1e4,5e4,2e3,1e3])
x = GenomicPositions( [ 200,100,150,50,500,20000], ["chr2","chr2","chr2","chrM","chrY","chr1"], chrinfo)
y = GenomicPositions( [ 220,100,50,50,420,10000], ["chr2","chr2","chr2","chrM","chrY","chr1"], chrinfo)
x in y
indexin(x,y)
slide!(x,20)
x in y
overlaps(x,y)

```


<a id='Ordering-1'></a>

## Ordering


Given the linearized representation of the genome, ordering operations are done with respect first to position within a chromosome and then to the order of chromosomes specified in the contained `GenomeInfo` object. Stand, for types that implement it, is not considered in ordering. This is different than how strand is handled in R's Bioconductor classes.


```julia
using DataFrames
chrinfo = GenomeInfo("GRCh38",["chr1","chr2","chr3","chrX","chrY","chrM"],Int64[3e5,2e5,1e4,5e4,2e3,1e3])
x = GenomicPositions( [ 200,100,150,50,500,20000], ["chr2","chr2","chr2","chrM","chrY","chr1"], chrinfo)
sortperm(x)
sort!(x,rev=true)
issorted(x)
y = sort(x)
issorted(y)
convert(DataFrame, y)
```


<a id='Intersection-/-overlap-operations-1'></a>

## Intersection / overlap operations


Currently we depend on Bio.jl and the `IntervalCollection` for overlap queries. We provide `convert` methods to make `IntervalCollection`s. These collections store the genome-scale positions put the genome string in the chromosome string field, resulting in a single tree. We add the index of each interval in the `AbstractGenomicVector` in the metadata slot of each `Interval` which can be used to relate the `IntervalCollection` back to our original object.

