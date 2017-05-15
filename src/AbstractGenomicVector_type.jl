##################################
### AbstractGenomicVector Type ###
##################################

abstract AbstractGenomicVector{T} <: AbstractVector{T}

"""AbstractGenomicVector API, requires

- findin
- indexin

"""
Returns an Logical Vector of positions in AbstratGenomicVector `x` that are in
AbstratGenomicVector `y` by exact match.
"""
function Base.in(x::AbstractGenomicVector, y::AbstractGenomicVector)
    inds = falses(length(x))
    inds[ findin(x,y) ] = true
    inds
end

"""
Like indexoverlap, but returns a Logical Vector indicating which of `x` have an overlap
in `y`.
"""
function hasoverlap(x::AbstractGenomicVector, y::AbstractGenomicVector)
  inds = falses(length(x))
  inds[ indexoverlap(x,y) ] = true
  inds
end

"""
Like intersect, but by overlap. (Should this be called subset_by_overlap?)
"""
overlap(x::AbstractGenomicVector, y::AbstractGenomicVector) = x[ indexoverlap(x,y) ]
