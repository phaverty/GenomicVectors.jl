__precompile__(true)

module GenomicVectors

# Re-implementation of the GenomicRanges type from Bioconductor's GenomicRanges package by H. Pages, P. Aboyoun and M.Lawrence

using RLEVectors
using AxisArrays
using GenomicFeatures
using DataFrames

# utils
export genopos

# types
import Base: ==, getindex, setindex!
import GenomicFeatures: Interval, IntervalCollection, Strand
export GenomeInfo, AbstractGenomicVector, GenomicPositions, GenomicRanges, GenomicDataFrame

# GenomeInfo
export chr_ends, chr_lengths, chr_offsets, chr_names, same_genome

# GenomicPositions and GenomicRanges
import Base: size, length, empty!, intersect, findin
import RLEVectors: starts, widths, ends, eachrange, each, disjoin
import GenomicFeatures: coverage
export genome, chr_info, chromosomes, chrpos, chrindex, genopos, slide!, slide
export strands, nearest, genostarts, genoends, starts, ends, widths
export _genostarts, _genoends, _strands
export findoverlaps, eachrange, remove_overlaps, select_overlaps, disjoin, gaps, coverage

# Delegations
import Base: similar, copy, unique, size, length, endof, issubset, vcat, union, intersect, setdiff, symdiff, append!, prepend!, resize!

# GenomicTable
export rowindex, table

### Includes
include("GenomeInfo_type.jl")
include("utils.jl")
include("AbstractGenomicVector_type.jl")
include("GenomicPositions_type.jl")
include("GenomicRanges_type.jl")
include("GenomicDataFrame_type.jl")
include("delegations.jl")
include("rcall.jl")
include("precompile.jl")

end # module
