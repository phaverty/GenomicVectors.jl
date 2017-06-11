##################################
### AbstractGenomicVector Type ###
##################################

## FIXME: Checking exact match possibly be efficient. Get Bio.Intervals to discriminate exact and overlap matching.

abstract AbstractGenomicVector{T} <: AbstractVector{T}

## Required API for an AbstractGenomicVector
## chr_info() must return a GenomeInfo object. This then allows all of the methods on GenomeInfo to work.
## _genostarts() must return a vector of start positions relative to the
## GenomeInfo. May be direct access, so do not modify values
## Similarly, _genoends() and _strands()

## Functions that work on the underlying GenomeInfo
for op in [:chr_names, :chr_lengths, :chr_ends, :chr_offsets, :genome]
    @eval $(op)(x::AbstractGenomicVector) = $(op)(chr_info(x))
end
"Return a Bool indicating if two objects represent positions on the identical genome."
same_genome(x::AbstractGenomicVector, y::AbstractGenomicVector) = chr_info(x) == chr_info(y)

## General purpose getters
"Get the starting nucleotide index for each range/position relative to the chromosome on which they lie."
starts(x::AbstractGenomicVector) = chrpos(_genostarts(x),chr_info(x))
"Get the ending nucleotide index for each range/position relative to the chromosome on which they lie."
ends(x::AbstractGenomicVector) = chrpos(_genoends(x),chr_info(x))
"Get the distance, between the start and end nucleotide of the range, 1 for positions."
widths(x::AbstractGenomicVector) = (_genoends(x) - _genostarts(x)) .+ 1
"Get the name of the chromosome for each range/position."
chromosomes(x::AbstractGenomicVector) = chromosomes(start(x),chr_info(x))
"Get the starting nucleotide index for each range/position in the linearized genome."
genostarts(x::AbstractGenomicVector) = copy(_genostarts(x))
"Get the ending nucleotide index for each range/position in the linearized genome."
genoends(x::AbstractGenomicVector) = copy(_genoends(x))
"Get the DNA strand for each range/position, pass by copy."
strands(x::AbstractGenomicVector) = copy(_strands(x))
"Return an iterator that returns tuples of genostart and genoend pairs."
each(x::AbstractGenomicVector) = zip(_genostarts(x),_genoends(x))

"""
The GenoPos Interface provides access to positional information in the linearized
genome or in chromosome coordinate (e.g. chr4:1000-1020).
"""
starts, ends, widths, chromosomes, genostarts, genoends, strands, each

"""
    findoverlaps(x::AbstractGenomicVector,y::AbstractGenomicVector)

Creates a `Bio.Intervals.IntersectIterator` from two `AbstractGenomicVectors`, much like the BioConductor
`findOverlaps`. `Bio.Intervals` calls this function `intersect`, but I would expect `intersect` to have the
same behavior as base, returning a subset copy of the first argument. `findoverlaps` is the kernel
of `findin`, `indexin` and `in`. 
"""
function findoverlaps(x::AbstractGenomicVector, y::AbstractGenomicVector)
    same_genome(x, y) || throw(ArgumentError("Both inputs must be from the same genome."))
    xit = convert(IntervalCollection,x)
    yit = convert(IntervalCollection,y)
    ol = intersect(xit,yit)
    ol
end

"""
    findin(x::AbstractGenomicVector,y::AbstractGenomicVector,exact::Bool=true)
    """
function Base.findin(x::AbstractGenomicVector, y::AbstractGenomicVector, exact::Bool=true)
    ol = findoverlaps(x,y)
    inds = Vector{Int64}(0)
    if exact
        for (el_a,el_b) in ol
            if exact_match(el_a,el_b)
                push!(inds,metadata(el_a))
            end
        end
    else
        for (el_a,el_b) in ol
            push!(inds,metadata(el_a))
        end
    end
    sort(unique(inds))
end

"""
    indexin(x::AbstractGenomicVector, y::AbstractGenomicVector, exact::Bool=true)
"""
function Base.indexin(x::AbstractGenomicVector, y::AbstractGenomicVector, exact::Bool=true)
    ol = findoverlaps(x,y)
    inds = zeros(Int64,length(x))
    if exact
        for (el_a,el_b) in ol
            if exact_match(el_a,el_b)
                inds[ metadata(el_a) ] = metadata(el_b)
            end
        end
    else
        for (el_a,el_b) in ol
            inds[ metadata(el_a) ] = metadata(el_b)
        end
    end
    inds
end
Base.in(x::AbstractGenomicVector, y::AbstractGenomicVector, exact::Bool=true) = indexin(x,y,exact) .!= 0
Base.intersect(x::AbstractGenomicVector, y::AbstractGenomicVector, exact::Bool=true) = x[ findin(x,y,exact) ]
Base.setdiff(x::AbstractGenomicVector, y::AbstractGenomicVector, exact::Bool=true) = x[!in(x,y,exact)]
exact_match(el_a::Interval, el_b::Interval) = first(el_a) == first(el_b) && last(el_a) == last(el_b)
"""
Matching functions in `GenomicVectors` can perform overlap matching, rather than exact
matching when given the extra argument `exact=false`. In either case, the genome strand
is never considered.
"""
findoverlaps, findin, indexin, in, intersect, setdiff

# The discussion in https://github.com/JuliaLang/Juleps/blob/master/Find.md#particular-cases is very relevant choosing names for these functions.
