abstract AbstractGenomicVector{T} <: AbstractVector{T}

## Required API for an AbstractGenomicVector
## chr_info() must return a GenomeInfo object.
##   This allows all of the methods on GenomeInfo to work.
## _genostarts() must return a vector of start positions relative to the GenomeInfo. May be direct access, so 

## Functions that work on the underlying GenomeInfo
## Documented by GenomeInfo
chr_info(x::AbstractGenomicVector) = x.chrinfo # Immutable, so just pass along
for op in [:chr_names, :chr_lengths, :chr_ends, :chr_offsets, :genome]
    @eval $(op)(x::AbstractGenomicVector) = $(op)(chr_info(x))
end

same_genome(x::AbstractGenomicVector, y::AbstractGenomicVector) = chr_info(x) == chr_info(y)
function same_genome(x::AbstractGenomicVector, y::Interval)
    seqname(y) != genome(x) && return false
    chrs = chromosomes( [leftposition(y),rightposition(y)], chr_info(x) )
    chrs[1] != chrs[2] && return false
    return true
end
same_genome(x::Interval, y::AbstractGenomicVector) = same_genome(y,x)

strands(x::AbstractGenomicVector) = RLEVector("+", length(x))

