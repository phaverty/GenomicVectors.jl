
<a id='The-GenomicVectors-Types-and-Methods-1'></a>

# The GenomicVectors Types and Methods


<a id='Index-1'></a>

## Index

- [`GenomicVectors.GenomeInfo`](api.md#GenomicVectors.GenomeInfo)
- [`GenomicVectors.GenomicPositions`](api.md#GenomicVectors.GenomicPositions)
- [`GenomicVectors.GenomicRanges`](api.md#GenomicVectors.GenomicRanges)
- [`GenomicVectors.chrpos`](api.md#GenomicVectors.chrpos)
- [`GenomicVectors.genopos`](api.md#GenomicVectors.genopos)
- [`GenomicVectors.nearest`](api.md#GenomicVectors.nearest)
- [`GenomicVectors.overlap`](api.md#GenomicVectors.overlap)
- [`GenomicVectors.overlapin`](api.md#GenomicVectors.overlapin)
- [`RLEVectors.ends`](api.md#RLEVectors.ends)
- [`RLEVectors.starts`](api.md#RLEVectors.starts)
- [`RLEVectors.widths`](api.md#RLEVectors.widths)


<a id='Types-1'></a>

## Types

<a id='GenomicVectors.GenomeInfo' href='#GenomicVectors.GenomeInfo'>#</a>
**`GenomicVectors.GenomeInfo`** &mdash; *Type*.



**GenomeInfo Type**

A GenomeInfo holds information about a genome including its name, chromosome names,  chromosome lengths and chromosome offsets into a concatenated, linear genome (genopos). Indexing returns the genopos end of the indexed chromosome.

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

<a id='GenomicVectors.GenomicPositions' href='#GenomicVectors.GenomicPositions'>#</a>
**`GenomicVectors.GenomicPositions`** &mdash; *Type*.



```
GenomicPositions(chrpos, chromosomes, genomeinfo)
GenomicPositions(genopos, genomeinfo)
```

Represents single-nucleotide positions in a genome.

This type uses its (immutable) `GenomeInfo` slot object to describe corresponding genome and     positions can be expressed relative to this concatenated, linearized genome or relative     to the chromosome containing a given position.

Sorting is by chromosome, as ordered by chrinfo,

By convention, all postions in a `GenomicPositions` are considered to be on the plus strand.

****

**Examples**

```julia
    genomeinfo = GenomeInfo("hg19",["chr1","chr2","chrX"],Int64[3e5,2e5,1e4])
    chrs = ["chr2","chr1","chr2","chrX"]
    pos = Int64[3e4,4.2e3,1.9e5,1e4]
    gpos = genopos(pos,chrs,chrinfo)
    x = GenomicPositions(pos,chrs,genomeinfo)
    y = GenomicPositions(gpos,genomeinfo)
    same_genome(x, y)
    sort!(y)
    convert(DataFrame, y)
```

<a id='GenomicVectors.GenomicRanges' href='#GenomicVectors.GenomicRanges'>#</a>
**`GenomicVectors.GenomicRanges`** &mdash; *Type*.



**`GenomicRanges`**

`GenomicRanges` represent closed ranges in a genome. This type uses its (immutable) `GenomeInfo` slot object to describe         corresponding genome and positions can be expressed relative to this concatenated, linearized genome or relative to the chromosome containing a given position.

**Examples**

```
chrinfo = GenomeInfo("hg19",["chr1","chr2","chrX"],Int64[3e5,2e5,1e4])
chrs = ["chr1","chr2","chr2","chrX"]
starts = [100, 200, 300, 400]
ends = [120, 240, 350, 455]
gr = GenomicRanges(chrs,starts,ends,chrinfo)
```

**Indexing**

Indexing a `GenomicRanges` with an array produces a new `GenomicRanges`. Indexing by a scalar produces a two-tuple of the start and end positions in genome location units.


<a id='Genome-Location-API-1'></a>

## Genome Location API


`GenomicVectors.jl` has ... All `AbstractGenomicVector`s implement the API for `GenomeInfo` for access to their genome descriptions.


<a id='Accessing-position-info-1'></a>

## Accessing position info


```
starts
widths
ends
genopos
chrpos
chromosomes
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


<a id='Overlap-function-1'></a>

### Overlap function


```
overlap
overlapin
overlapindex
nearest
```

