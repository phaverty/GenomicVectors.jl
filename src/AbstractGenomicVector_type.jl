##################################
### AbstractGenomicVector Type ###
##################################

## FIXME: Checking exact match possibly be efficient. Get Bio.Intervals to discriminate exact and overlap matching.

abstract AbstractGenomicVector{T} <: AbstractVector{T}

"""AbstractGenomicVector API, requires

- findin
- indexin

"""
function Base.findin(x::AbstractGenomicVector, y::AbstractGenomicVector)
    same_genome(x, y) || throw(ArgumentError("Both inputs must be from the same genome."))
    xit = convert(IntervalCollection,x)
    yit = convert(IntervalCollection,y)
    ol = intersect(xit,yit)
    inds = Vector{Int64}(0)
    for (el_a,el_b) in ol
        if first(el_a) == first(el_b) && last(el_a) == last(el_b) && strand(el_a) == strand(el_b)
            push!(inds,metadata(el_a))
        end
    end
    sort(unique(inds))
end

function Base.indexin(x::AbstractGenomicVector, y::AbstractGenomicVector)
    same_genome(x, y) || throw(ArgumentError("Both inputs must be from the same genome."))
    xit = convert(IntervalCollection,x)
    yit = convert(IntervalCollection,y)
    ol = intersect(xit,yit)
    inds = zeros(Int64,3)
    for (i,(el_a,el_b)) in enumerate(ol)
        if first(el_a) == first(el_b) && last(el_a) == last(el_b) && strand(el_a) == strand(el_b)
            inds[i] = metadata(el_b)
        else
            inds[i] = 0
        end
    end
    inds
end

Base.in(x::AbstractGenomicVector, y::AbstractGenomicVector) = indexin(x,y) .!= 0
Base.intersect(x::AbstractGenomicVector, y::AbstractGenomicVector) = x[ findin(x,y) ]
Base.setdiff(x::AbstractGenomicVector, y::AbstractGenomicVector) = x[!in(x,y)]
