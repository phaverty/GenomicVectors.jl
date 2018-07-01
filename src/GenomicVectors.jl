__precompile__(true)

module GenomicVectors

# Re-implementation of the GenomicRanges type from Bioconductor's GenomicRanges package by H. Pages, P. Aboyoun and M.Lawrence

using RLEVectors
using AxisArrays
using GenomicFeatures
using DataFrames
using BioAlignments
#using RCall

# utils
export genopos

# types
import Base: ==, getindex, setindex!
import GenomicFeatures: Interval, IntervalCollection, Strand, strand, leftposition, rightposition
export GenomeInfo, AbstractGenomicVector, GenomicPositions, GenomicRanges
export GenomicDataFrame, GenomicVectorIterator

# GenomeInfo
export chr_ends, chr_lengths, chr_offsets, chr_names, same_genome

# GenomicPositions and GenomicRanges
import Base: size, length, empty!, intersect, findin
import RLEVectors: starts, widths, ends, eachrange, each, disjoin
import GenomicFeatures: coverage
export genome, chr_info, chromosomes, chrpos, chrindex, genopos, slide!, slide
export strands, nearest, genostarts, genoends, starts, ends, widths
export _genostarts, _genoends, _strands
export findoverlaps, overlap_table, eachrange
#export remove_overlaps, select_overlaps
export disjoin, gaps, coverage, collapse

# Delegations
import Base: similar, copy, unique, size, length, endof, issubset, vcat, union, intersect, setdiff, symdiff, append!, prepend!, resize!

# GenomicTable
export rowindex, table

# RCall
#import RCall: sexp, rcopy, RClass, rcopytype, @R_str, S4Sxp

# BAM
export strand

### Includes
include("GenomeInfo_type.jl")
include("utils.jl")
include("AbstractGenomicVector_type.jl")
include("GenomicPositions_type.jl")
include("GenomicRanges_type.jl")
include("GenomicDataFrame_type.jl")
include("delegations.jl")
#include("rcall.jl")
include("bam.jl")
include("precompile.jl")
_precompile_()

end # module
