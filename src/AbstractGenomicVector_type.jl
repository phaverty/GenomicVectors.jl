##################################
### AbstractGenomicVector Type ###
##################################

"""
An AbstractGenomicVector is a Vector that describes positions or ranges in a
single genome, in an arbitrary order. An AbstractGenomicVector must implement
the GenomeInfo and GenoPos Interfaces. Sorting is by chromosome then by
nucleotide position.
"""
abstract AbstractGenomicVector{T} <: AbstractVector{T}

Base.sort(x::AbstractGenomicVector; rev::Bool=false) = sort!(copy(x))
Base.issorted(x::AbstractGenomicVector; rev::Bool=false) = issorted( eachrange(x), rev=rev )
Base.sortperm(x::AbstractGenomicVector; rev=false) = sortperm( collect(eachrange(x)), rev=rev ) # No method for iterator
_exact_overlap(el_a::Interval, el_b::Interval) = first(el_a) == first(el_b) && last(el_a) == last(el_b)
slide(x::AbstractGenomicVector, value::Integer) = slide!( copy(x), value )

function Base.show(io::IO, x::AbstractGenomicVector)
    Base.show_vector(io,convert(Vector{String},x),"[", "]")
end    

function Base.show(io::IO, ::MIME"text/plain", x::AbstractGenomicVector)
    t = typeof(x)::DataType
    show(io, t)
    show(io, convert(GenomicTable, x))
end

if VERSION >= v"0.6.0"
    Base.IndexStyle(::Type{<:AbstractGenomicVector}) = IndexLinear()
else
    Base.linearindexing{T<:AbstractGenomicVector}(::Type{T}) = Base.LinearFast()
end

function findoverlaps(x::AbstractGenomicVector, y::AbstractGenomicVector, exact::Bool=false)
    same_genome(x, y) || throw(ArgumentError("Both inputs must be from the same genome."))
    xit = convert(IntervalCollection,x)
    yit = convert(IntervalCollection,y)
    if exact
        ol = eachoverlap(xit, yit, filter=_exact_overlap)
    else
        ol = eachoverlap(xit,yit)
    end
    ol
end

function Base.findin(x::AbstractGenomicVector, y::AbstractGenomicVector, exact::Bool=true)
    ol = findoverlaps(x,y)
    inds = Vector{Int64}(0)
    for (el_a,el_b) in ol
        push!(inds,metadata(el_a))
    end
    sort(unique(inds))
end

function Base.indexin(x::AbstractGenomicVector, y::AbstractGenomicVector, exact::Bool=true)
    ol = findoverlaps(x,y)
    inds = zeros(Int64,length(x))
    for (el_a,el_b) in ol
        inds[ metadata(el_a) ] = metadata(el_b)
    end
    inds
end
Base.in(x::AbstractGenomicVector, y::AbstractGenomicVector, exact::Bool=true) = indexin(x,y,exact) .!= 0
Base.intersect(x::AbstractGenomicVector, y::AbstractGenomicVector, exact::Bool=true) = x[ findin(x,y,exact) ]
Base.setdiff(x::AbstractGenomicVector, y::AbstractGenomicVector, exact::Bool=true) = x[!in(x,y,exact)]

"""
# Searching GenomicVectors

Matching functions in `GenomicVectors` can perform overlap matching, rather than exact
matching when given the extra argument `exact=false`. In either case, the genome strand
is never considered.

    findoverlaps(x::AbstractGenomicVector,y::AbstractGenomicVector)

Creates a `Bio.Intervals.IntersectIterator` from two `AbstractGenomicVectors`, much like the BioConductor
`findOverlaps`. `Bio.Intervals` calls this function `intersect`, but I would expect `intersect` to have the
same behavior as base, returning a subset copy of the first argument. `findoverlaps` is the kernel
of the other search/set operations.

    findin(x::AbstractGenomicVector,y::AbstractGenomicVector,exact::Bool=true)

    indexin(x::AbstractGenomicVector, y::AbstractGenomicVector, exact::Bool=true)

    in(x::AbstractGenomicVector, y::AbstractGenomicVector, exact::Bool=true)
The `AbstractGenomicVector` method on `in` is vectorized and returns a `BitArray`
that is `true` for each element of `x` that is in the set `y`.
    
    intersect(x::AbstractGenomicVector, y::AbstractGenomicVector, exact::Bool=true)

    setdiff(x::AbstractGenomicVector, y::AbstractGenomicVector, exact::Bool=true)
    
"""
findoverlaps, findin, indexin, in, intersect, setdiff

# The discussion in https://github.com/JuliaLang/Juleps/blob/master/Find.md#particular-cases is very relevant choosing names for these functions.
