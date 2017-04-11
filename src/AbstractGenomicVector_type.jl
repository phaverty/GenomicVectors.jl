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

## General purpose getters
function chromosomes(x::AbstractGenomicVector)
    names = chr_names(x)
    ends = chr_ends(x.chrinfo)
    offsets = chr_offsets(x.chrinfo)
    nchr = length(names)
    res = similar(names, length(x))
    r = 1
    i = 1
    @inbounds for pos in _genostarts(x)
        if pos > ends[r] || pos <= offsets[r]
            r = searchsortedfirst(ends, pos, 1, nchr, Base.Forward)
        end
        res[i] = names[r]
        i = i + 1
    end
    res
end
strands(x::AbstractGenomicVector) = RLEVector("+", length(x))

