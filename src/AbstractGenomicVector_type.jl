##################################
### AbstractGenomicVector Type ###
##################################

## FIXME: Checking exact match possibly be efficient. Get Bio.Intervals to discriminate exact and overlap matching.

"""
An AbstractGenomicVector is a Vector that descrbes positions or ranges in a
single genome, in an arbitrary order. An AbstractGenomicVector must implement
the GenomeInfo and GenoPos Interfaces.
"""
abstract AbstractGenomicVector{T} <: AbstractVector{T}

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
            if _exact_match(el_a,el_b)
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
            if _exact_match(el_a,el_b)
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
_exact_match(el_a::Interval, el_b::Interval) = first(el_a) == first(el_b) && last(el_a) == last(el_b)

"""
Matching functions in `GenomicVectors` can perform overlap matching, rather than exact
matching when given the extra argument `exact=false`. In either case, the genome strand
is never considered.
"""
findoverlaps, findin, indexin, in, intersect, setdiff

# The discussion in https://github.com/JuliaLang/Juleps/blob/master/Find.md#particular-cases is very relevant choosing names for these functions.
