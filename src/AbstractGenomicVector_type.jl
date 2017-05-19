##################################
### AbstractGenomicVector Type ###
##################################

## FIXME: Checking exact match possibly be efficient. Get Bio.Intervals to discriminate exact and overlap matching.

abstract AbstractGenomicVector{T} <: AbstractVector{T}

"""AbstractGenomicVector API

    findoverlaps(agv1, agv2)
Creates a `Bio.Intervals.IntersectIterator` from two `AbstractGenomicVectors`, much like the BioConductor
`findOverlaps`. `Bio.Intervals` calls this function `intersect`, but I would expect `intersect` to have the
    same behavior as base, returning a subset copy of the first argument. `findoverlaps` is the kernel
    of `findin`, `indexin` and `in`.

The discussion in https://github.com/JuliaLang/Juleps/blob/master/Find.md#particular-cases is very relevant
    to the naming issue.
"""
function findoverlaps(x::AbstractGenomicVector, y::AbstractGenomicVector)
    same_genome(x, y) || throw(ArgumentError("Both inputs must be from the same genome."))
    xit = convert(IntervalCollection,x)
    yit = convert(IntervalCollection,y)
    ol = intersect(xit,yit)
    ol
end

function Base.findin(x::AbstractGenomicVector, y::AbstractGenomicVector, exact::Bool=true)
    ol = findoverlaps(x,y)
    inds = Vector{Int64}(0)
    if exact
        for (el_a,el_b) in ol
            if first(el_a) == first(el_b) && last(el_a) == last(el_b) && strand(el_a) == strand(el_b)
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

function Base.indexin(x::AbstractGenomicVector, y::AbstractGenomicVector, exact::Bool=true)
    ol = findoverlaps(x,y)
    inds = zeros(Int64,length(x))
    if exact
        for (i,(el_a,el_b)) in enumerate(ol)
            if first(el_a) == first(el_b) && last(el_a) == last(el_b) && strand(el_a) == strand(el_b)
                inds[ metadata(el_a) ] = metadata(el_b)
            end
        end
    else
        for (i,(el_a,el_b)) in enumerate(ol)
                inds[ metadata(el_a) ] = metadata(el_b)
        end
    end
    inds
end
Base.in(x::AbstractGenomicVector, y::AbstractGenomicVector, exact::Bool=true) = indexin(x,y,exact) .!= 0
Base.intersect(x::AbstractGenomicVector, y::AbstractGenomicVector, exact::Bool=true) = x[ findin(x,y,exact) ]
Base.setdiff(x::AbstractGenomicVector, y::AbstractGenomicVector, exact::Bool=true) = x[!in(x,y,exact)]
