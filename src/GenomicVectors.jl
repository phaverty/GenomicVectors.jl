VERSION >= v"0.5.0" && __precompile__(true)

module GenomicVectors

# Re-implementation of the GenomicRanges type from Bioconductor's GenomicRanges package by H. Pages, P. Aboyoun and M.Lawrence

using RLEVectors
using AxisArrays
using Bio.Intervals
using DataFrames

# utils
export genopos

# types
import Base: ==, getindex, setindex!
import Bio.Intervals: StringField, Interval, IntervalCollection, Strand
export GenomeInfo, AbstractGenomicVector, GenomicPositions, GenomicRanges

# GenomeInfo
export chr_ends, chr_lengths, chr_offsets, chr_names, same_genome

# GenomicPositions and GenomicRanges
import Base: size, length, empty!, intersect, findin
import RLEVectors: starts, widths, ends, each
export genome, chr_info, chromosomes, chrpos, genopos, slide!, slide, in, overlaps, starts, ends, widths, strands, indexin, overlapindex, nearest, genostarts, genoends
export getindex, setindex!, overlapin, overlap, in, each, hasoverlap

# Delegations
import Base: similar, copy, unique, size, length, endof, issubset, vcat, union, intersect, setdiff, symdiff, append!, prepend!
export similar, copy, unique, size, length, endof, issubset, vcat, union, intersect, setdiff, symdiff, append!, prepend!

### Includes
include("GenomeInfo_type.jl")
include("AbstractGenomicVector_type.jl")
include("utils.jl")
include("GenomicPositions_type.jl")
include("GenomicRanges_type.jl")
include("delegations.jl")

if VERSION >= v"0.5.0"
    include("precompile.jl")
    _precompile_()
end


end # module
