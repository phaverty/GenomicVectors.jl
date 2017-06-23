VERSION >= v"0.5.0" && __precompile__(true)

module GenomicVectors

# Re-implementation of the GenomicRanges type from Bioconductor's GenomicRanges package by H. Pages, P. Aboyoun and M.Lawrence

using RLEVectors
using AxisArrays
using GenomicFeatures
using DataTables

# utils
export genopos

# types
import Base: ==, getindex, setindex!
import GenomicFeatures: Interval, IntervalCollection, Strand
export GenomeInfo, AbstractGenomicVector, GenomicPositions, GenomicRanges, GenomicTable

# GenomeInfo
export chr_ends, chr_lengths, chr_offsets, chr_names, same_genome

# GenomicPositions and GenomicRanges
import Base: size, length, empty!, intersect, findin
import RLEVectors: starts, widths, ends, eachrange
export genome, chr_info, chromosomes, chrpos, chrindex, genopos, slide!, slide, strands, nearest, genostarts, genoends, starts, ends, widths, findoverlaps

# Delegations
import Base: similar, copy, unique, size, length, endof, issubset, vcat, union, intersect, setdiff, symdiff, append!, prepend!

# GenomicTable
export rowindex, table

### Includes
include("GenomeInfo_type.jl")
include("utils.jl")
include("AbstractGenomicVector_type.jl")
include("GenomicPositions_type.jl")
include("GenomicRanges_type.jl")
include("GenomicTable_type.jl")
include("delegations.jl")

if VERSION >= v"0.5.0"
    include("precompile.jl")
    _precompile_()
end


end # module
