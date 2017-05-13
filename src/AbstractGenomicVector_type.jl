### We will try to provide much of the functionality of the types in this
### package via the GenomeInfo Interface, the GenoPos Interface and methods
### for the AbstractGenomicVector type.

### The GenomeInfo Interface provides
### access to the genome (name, chromosome lengths, etc.). This Interface requires
### the following methods:
### - chr_info
### which returns a GenomeInfo object, which then provides these methods:
### - chr_names
### - chr_lengths
### - chr_ends
### - chr_offsets
### - genome
### - same_genome

### The GenoPos Interface provides access to positional information in the linearized
### genome or in chromome coordinate (e.g. chr4:1000-1020).
abstract AbstractGenomicVector{T} <: AbstractVector{T}

for op in [:chr_names, :chr_lengths, :chr_ends, :chr_offsets, :genome]
    @eval $(op)(x) = $(op)(chr_info(x))
end

same_genome(x, y) = chr_info(x) == chr_info(y)
function same_genome(x, y::Interval)
    seqname(y) != genome(x) && return false
    chrs = chromosomes( [leftposition(y),rightposition(y)], chr_info(x) )
    chrs[1] != chrs[2] && return false
    return true
end
same_genome(x::Interval, y) = same_genome(y,x)

@doc (@doc AbstractGenomicVector) GenomeInfo, chr_info, chr_names, chr_lengths, chr_ends, chr_offsets, genome, same_genome
