######################
### GenomicRanges Type
######################

"""
# `GenomicRanges`
`GenomicRanges` represent closed ranges in a genome. This type uses its (immutable) `GenomeInfo` slot object to describe
        corresponding genome and positions can be expressed relative to this concatenated, linearized genome or relative to the chromosome containing a given position.

## Examples
    chrinfo = GenomeInfo("hg19",["chr1","chr2","chrX"],Int64[3e5,2e5,1e4])
    chrs = ["chr1","chr2","chr2","chrX"]
    starts = [100, 200, 300, 400]
    ends = [120, 240, 350, 455]
    gr = GenomicRanges(chrs,starts,ends,chrinfo)

## Indexing
Indexing a `GenomicRanges` with an array produces a new `GenomicRanges`. Indexing by a scalar produces a two-tuple of the start and end positions in genome location units.

"""
type GenomicRanges{T1 <: Integer} <: AbstractGenomicVector{T1}
    starts::Vector{T1}
    ends::Vector{T1}
    strands::Vector{Strand}
    chrinfo::GenomeInfo{T1}
    function GenomicRanges(starts, ends, strands, chrinfo)
        length(starts) != length(ends) && throw(ArgumentError("starts and ends must be of the same length."))
        if strands == nothing
            strands = Vector{Strand}(length(starts))
            strands[:] = STRAND_NA
        else
            length(starts) != length(strands) && throw(ArgumentError("starts, ends and stands must be of the same length."))
        end
        new(starts, ends, strands, chrinfo)
    end
end
## Create with specified strands
GenomicRanges{T1 <: Integer}(chrs::Vector{String}, starts::Vector{T1}, ends::Vector{T1}, strands::Vector{Char}, chrinfo::GenomeInfo{T1}) = GenomicRanges{T1}(genopos(starts,chrs,chrinfo), genopos(ends,chrs,chrinfo),strands,chrinfo)
GenomicRanges{T1 <: Integer}(chrs::Vector{String}, starts::Vector{T1}, ends::Vector{T1}, strands::Vector{Strand}, chrinfo::GenomeInfo{T1}) = GenomicRanges{T1}(genopos(starts,chrs,chrinfo), genopos(ends,chrs,chrinfo),strands,chrinfo)
GenomicRanges{T1 <: Integer}(genostarts::Vector{T1}, genoends::Vector{T1}, strands::Vector{Char}, chrinfo::GenomeInfo{T1}) = GenomicRanges{T1}(genostarts,genoends,strands,chrinfo)
GenomicRanges{T1 <: Integer}(genostarts::Vector{T1}, genoends::Vector{T1}, strands::Vector{Strand}, chrinfo::GenomeInfo{T1}) = GenomicRanges{T1}(genostarts,genoends,strands,chrinfo)
## Create with default strands
GenomicRanges{T1 <: Integer}(chrs::Vector{String}, starts::Vector{T1}, ends::Vector{T1}, chrinfo::GenomeInfo{T1}) = GenomicRanges{T1}(genopos(starts,chrs,chrinfo), genopos(ends,chrs,chrinfo), nothing, chrinfo)
GenomicRanges{T1 <: Integer}(genostarts::Vector{T1}, genoends::Vector{T1}, chrinfo::GenomeInfo{T1}) = GenomicRanges{T1}(genostarts,genoends,nothing,chrinfo)

## Getters
_genostarts(x::GenomicRanges) = x.starts # Pass by reference for internal use
_genoends(x::GenomicRanges) = x.ends # Pass by reference for internal use
_strands(x::GenomicRanges) = x.strands # Pass by reference for internal use

genostarts(x::GenomicRanges) = copy(x.starts)
genoends(x::GenomicRanges) = copy(x.ends)
strands(x::GenomicRanges) = copy(x.strands)

starts(x::GenomicRanges) = chrpos(x.starts,chr_info(x))
ends(x::GenomicRanges) = chrpos(x.ends,chr_info(x))
widths(x::GenomicRanges) = (x.ends - x.starts) .+ 1
chromosomes(x::GenomicRanges) = chromosomes(_genostarts(x),chr_info(x))

## Indexing
each(x::GenomicRanges) = zip(x.starts,x.ends)
Base.getindex(x::GenomicRanges, i::Int) = (x.starts[i],x.ends[i],x.strands[i])

function Base.getindex(x::GenomicRanges, i::AbstractArray)
    GenomicRanges( x.starts[i], x.ends[i], x.strands[i], chr_info(x) )
end

function Base.setindex!(x::GenomicRanges, value, i)
    (min,max) = extrema(i)
    if min < 1 || max > x.chrinfo.chr_ends[end]
        error("Incoming genopos is outside the bounds of the genome.")
    end
    x.genopos[i] = value
    return(x)
end

function Base.vcat(x::GenomicRanges,y::GenomicRanges)
    same_genome(x, y) || throw(ArgumentError("Both GenomicPositions must be from the same genome."))
    GenomicRanges(vcat(x.starts,y.starts),vcat(x.ends,y.ends),vcat(x.strands,y.strands),chr_info(x))
end

## Show
function Base.show(io::IO, ::MIME"text/plain", x::GenomicRanges)
    t = typeof(x)::DataType
    show(io, t)
    write(io, "\nGenome Metadata:\n  ")
    show(io, x.chrinfo)
    write(io, "\nChromosomes:\n ")
    Base.show_vector(io, chromosomes(x),"[","]")
    write(io, "\nChromosome Start Positions:\n ")
    Base.show_vector(io, starts(x),"[","]")
    write(io, "\nChromosome End Positions:\n ")
    Base.show_vector(io, ends(x),"[","]")
    write(io, "\nChromosome Strand:\n ")
    Base.show_vector(io, strands(x),"[","]")
end

## Conversions
function Base.convert(::Type{DataFrame}, x::GenomicRanges)
    n = length(x)
    chrs = chr_names(x)
    n_chrs = length(chrs)
    c = similar(chrs, n)
    s = similar(x.starts, n)
    e = similar(x.ends, n)
    ends = chr_ends(x.chrinfo)
    offsets = chr_offsets(x.chrinfo)
    i = 1
    for (spos,epos) in zip(x.starts,x.ends)
        ind = searchsortedfirst(ends, spos, one(Int64), n_chrs, Base.Forward)
        c[i] = chrs[ ind ]
        o = offsets[ ind ]
        s[i] = spos - o
        e[i] = epos - o
        i = i + 1
    end
    return( DataFrame( Chromosome=c, Start=s, End=e, Strand=copy(strands(x)) ) )
end

function Base.convert(::Type{Vector{String}}, x::GenomicRanges)
    df = convert(DataFrame,x)
    String[ string(c, ":", s, "-", e) for (c,s,e) in zip(df[:Chromosome], df[:Start], df[:End]) ]
end

Base.convert(::Type{Vector}, x::GenomicRanges) = [ (s,e) for (s,e) in x ]

Base.convert(::Type{GenomicPositions}, x::GenomicRanges) = GenomicPositions(starts(x), chr_info(x))

"""
Conversion of GenomicRanges to IntervalCollection adds index as metadata in order to recover order later.
"""
function Base.convert(::Type{IntervalCollection}, x::GenomicRanges)
    g = genome(x)
    IntervalCollection( sort([Interval(g,b,e,s,i) for (i,(b,e,s)) in enumerate(x)]) )
end

## Altering Positions
function slide!(gr::GenomicRanges, x::Integer)
    offsets = chr_offsets(gr)
    ends = chr_ends(gr)
    n_chrs = length(ends)
    chr_ind = 1
    i = 1
    for (s,e) in gr
        if e > ends[chr_ind] || s <= offsets[chr_ind] # Find new chr
            chr_ind = searchsortedfirst(ends, s, one(Int64), n_chrs, Base.Forward)
        end
        news = s + x
        newe = e + x
        if newe > ends[chr_ind] || news <= offsets[chr_ind]
            throw(ArgumentError("Genomic position ($s,$e) falls outside the bounds of chromosome $(chr_names(gr)[chr_ind]) when shifted by $x."))
        end
        gr.starts[i] = news
        gr.ends[i] = newe
        i = i + 1
    end
    gr
end

slide(gr::GenomicRanges, x::Integer) = slide!( copy(gr), x )

function Base.empty!(x::GenomicRanges)
    empty!(x.starts)
    empty!(x.ends)
    empty!(x.strands)
    x
end

## Sorting
function Base.sort!(x::GenomicRanges; rev::Bool=false)
    mat = sortrows( hcat(x.starts, x.ends), rev=rev )
    x.starts = mat[:,1]
    x.ends = mat[:,2]
    x
end

function Base.sort(x::GenomicRanges; rev::Bool=false)
    mat = sortrows( hcat(x.starts, x.ends), rev=rev )
    GenomicRanges( mat[:,1], mat[:,2], chr_info(x) )
end

function Base.issorted(x::GenomicRanges; rev::Bool=false)
    length(x) == 0 && return true
    (prev_s,prev_e) = x[1]
    if rev
        for (s,e) in x
            if s > prev_s || (s == prev_s && e > prev_e)
                return false
            end
        end
    else
        for (s,e) in x
            if s < prev_s || (s == prev_s && e < prev_e)
                return false
            end
        end
    end
    true
end

function Base.sortperm(x::GenomicRanges; rev=false)
    sortperm( convert(Vector,gr), rev=rev )
end

## Querying Positions
# Note that the standard set operations require exact matches and
# a separate set of functions work on overlaps
# N.B. intersect(tree1,tree2) returns an iterator of overlapping pairs. This is like
# The twocol matrix of match pairs in R's IRanges::findOverlaps()

# Identical matches (set ops)
function Base.findin(x::GenomicRanges, y::GenomicRanges)
    same_genome(x, y) || throw(ArgumentError("Both GenomicPositions must be from the same genome."))
    ## FIXME: This can't possibly be efficient. Get Bio.Intervals to discriminate exact and overlap matching.
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

function Base.in(x::GenomicRanges, y::GenomicRanges)
    inds = falses(length(x))
    inds[ findin(x,y) ] = true
    inds
end

Base.setdiff(x::GenomicRanges, y::GenomicRanges) = x[ !in(x,y) ]
Base.intersect(x::GenomicRanges, y::GenomicRanges) = x[ findin(x,y) ]

function Base.union(x::GenomicRanges, y::GenomicRanges)
    ic = convert(IntervalCollection,vcat(x,y))
    ic = unique(ic)
    inds = [ metadata(el) for el in ic ]
    gr = GenomicRanges( [ first(el) for el in ic ], [ last(el) for el in ic ], [ strand(el) for el in ic ], chr_info(x) )
    gr = gr[inds]
    gr
end

# Overlap ops
"""
Finds sorted, unique indexes of `x` that are in `y`. In other words, it i like
findin, but by range overlap.
"""
function overlapin(x::GenomicRanges, y::GenomicRanges)
    same_genome(x, y) || throw(ArgumentError("Both GenomicPositions must be from the same genome."))
    xit = convert(IntervalCollection,x)
    yit = convert(IntervalCollection,y)
    ol = intersect(xit,yit)
    sort( unique( [ metadata(el_a) for (el_a,el_b) in ol ] ) )
end

"""
Like overlapin, but returns a BitVector indicating which of `x` have an overlapin
in `y`.
"""
function hasoverlap(x::GenomicRanges, y::GenomicRanges)
  inds = falses(length(x))
  inds[ overlapin(x,y) ] = true
  inds
end

"""
Like intersect, but by overlap.
"""
overlap(x::GenomicRanges, y::GenomicRanges) = x[ overlapin(x,y) ]
