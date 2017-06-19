#########################
### GenomicPositions Type
#########################

"""
# GenomicPositions Type

Represents single-nucleotide positions in a genome.

This type uses its (immutable) `GenomeInfo` slot object to describe corresponding
genome. Therefore, positions can be expressed relative to this concatenated,
linearized genome or relative to the chromosome containing a given position.

By convention, all postions in a `GenomicPositions` are considered to be on the plus strand.

# Examples
```julia
    genomeinfo = GenomeInfo("hg19",["chr1","chr2","chrX"],Int64[3e5,2e5,1e4])
    chrs = ["chr2","chr1","chr2","chrX"]
    pos = Int64[3e4,4.2e3,1.9e5,1e4]
    gpos = genopos(pos,chrs,chrinfo)
    x = GenomicPositions(pos,chrs,genomeinfo)
    y = GenomicPositions(gpos,genomeinfo)
    same_genome(x, y)
    convert(DataTable, y)
```
"""
immutable GenomicPositions{T1 <: Integer} <: AbstractGenomicVector{T1}
    genopos::Vector{T1}
    chrinfo::GenomeInfo{T1}
    function GenomicPositions(genopos, chrinfo)
        new(genopos,chrinfo)
    end
end
GenomicPositions{T1 <: Integer}(genopos::Vector{T1}, chrinfo::GenomeInfo{T1}) = GenomicPositions{T1}(genopos, chrinfo)
GenomicPositions{T1 <: Integer}(pos::Vector{T1}, chromosomes::Vector{String}, chrinfo::GenomeInfo{T1}) = GenomicPositions{T1}(genopos(pos, chromosomes, chrinfo),chrinfo)

## GenomeInfo Interface
chr_info(x::GenomicPositions) = x.chrinfo

## GenoPos Interface
_genostarts(x::GenomicPositions) = x.genopos # Pass by reference for internal use
_genoends(x::GenomicPositions) = x.genopos # Pass by reference for internal use
_strands(x::GenomicPositions) = RLEVector("+", length(x))
strands(x::GenomicPositions) = _strands(x)
widths(x::GenomicPositions) = RLEVector(1, length(x))

## Indexing
Base.getindex(x::GenomicPositions, i::Int) = x.genopos[i]
Base.getindex(x::GenomicPositions, i::AbstractVector) = GenomicPositions( x.genopos[i], x.chrinfo )

function Base.setindex!(x::GenomicPositions, value, i)
    (min,max) = extrema(i)
    if min < 1 || max > x.chrinfo.chr_ends[end]
        throw(BoundsError("Incoming genopos is outside the bounds of the genome."))
    end
    x.genopos[i] = value
    return(x)
end

## Show
function Base.show(io::IO, x::GenomicPositions)
    t = typeof(x)::DataType
    show(io, t)
    write(io, "\nGenome Metadata:\n  ")
    show(io, x.chrinfo)
    write(io, "\nChromosomes:\n  ")
    Base.show_vector(io, chromosomes(x),"[","]")
    write(io, "\nChromosome Positions:\n  ")
    Base.show_vector(io, starts(x),"[","]")
end

## Conversions
function Base.convert(::Type{DataTable}, x::GenomicPositions)
    chrs = chr_names(x)
    n = length(x)
    c_res = similar(chrs, n)
    p_res = similar(_genostarts(x), n)
    ends = chr_ends(x)
    offsets = chr_offsets(x)
    i = r = 1
    e = ends[r]
    o = offsets[r]
    c = chrs[r]
    @inbounds for g in _genostarts(x)
        if g > e || g <= o
            r = 1
            while g > ends[r]
                r = r + 1
            end
            e = ends[r]
            o = offsets[r]
            c = chrs[r]
        end
        c_res[i] = c
        p_res[i] = g - o
        i = i + 1
    end
    return( DataTable( [c_res,p_res], [:Chromosome, :Position] ) )
end

function Base.convert(::Type{Vector{String}}, x::GenomicPositions)
    df = convert(DataTable,x)
    eltype(chr_names(x))[ string(chr, ":", pos, "-", pos) for (chr,pos) in zip(df[:Chromosome], df[:Position]) ]
end
Base.convert(::Type{Vector}, x::GenomicPositions) = genostarts(x)

"""
Conversion of GenomicPositions to IntervalCollection adds index as metadata in order to recover order later.
"""
function Base.convert(::Type{IntervalCollection}, x::GenomicPositions)
    g = genome(x)
    IntervalCollection( sort([Interval(g,s,s,STRAND_POS,i) for (i,s) in enumerate(x)]) )
end

## Altering Positions
function slide!(gpos::GenomicPositions, x::Integer)
    offsets = chr_offsets(gpos)
    ends = chr_ends(gpos)
    n_chrs = length(ends)
    chr_ind = 1
    i = 1
    for g in _genostarts(gpos)
        if g > ends[chr_ind] || g <= offsets[chr_ind] # Find new chr
            chr_ind = searchsortedfirst(ends, g, one(Int64), n_chrs, Base.Forward)
        end
        newg = g + x
        if newg > ends[chr_ind] || newg <= offsets[chr_ind]
            throw(ArgumentError("Genomic position $g falls outside the bounds of chromosome $(chr_names(gpos)[chr_ind]) when shifted by $x."))
        end
        gpos.genopos[i] = newg
        i = i + 1
    end
    gpos
end

slide(gr::GenomicPositions, x::Integer) = slide!( copy(gr), x )

Base.empty!(x::GenomicPositions) = empty!(x.genopos)

## Sorting
Base.sort(x::GenomicPositions; rev::Bool=false) = GenomicPositions( sort(genostarts(x), rev=rev), chr_info(x) )
Base.sort!(x::GenomicPositions; rev::Bool=false) = sort!(x.genopos, rev=rev)
Base.issorted(x::GenomicPositions; rev=false) = issorted(genostarts(x), rev=rev)
Base.sortperm(x::GenomicPositions; rev=false) = sortperm(genostarts(x), rev=rev)
Base.in(query::GenomicPositions, target::GenomicPositions, exact::Bool=true) = [in(v,target) for v in query]

"""
    function nearest(query::GenomicPositions, target::GenomicPositions)
For each `query` finds index in `target` that is nearest on the same chromosome.
If no match on the same chromosome exists, the index will be 0.
"""
function nearest{T}(query::GenomicPositions{T}, target::GenomicPositions{T})
    same_genome(query, target) || throw(ArgumentError("query and target must be from the same genome."))
    target_gpos = target.genopos
    target_chrs = chromosomes(target)
    nquery = length(query)
    res = Vector{Int64}(nquery)
    i = 1
    for (qpos, qchr) in zip(genostarts(query), chromosomes(query))
        temp_min = typemax(Int64)
        temp_min_index = 0
        j = 1
        for (tpos, tchr) in zip(target_gpos, target_chrs)
            new_min = abs(tpos - qpos)
            if new_min < temp_min && qchr == tchr
                temp_min = new_min
                temp_min_index = j
            end
            j = j + 1
        end
        res[i] = temp_min_index
        i = i + 1
    end
    res
end
