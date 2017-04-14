
"""
# GenomeInfo Type
A GenomeInfo holds information about a genome including its name, chromosome names,
 chromosome lengths and chromosome offsets into a concatenated, linear genome (genopos).
Indexing returns the genopos end of the indexed chromosome.

## Examples
```julia
chrinfo = GenomeInfo("hg19",["chr1","chr2","chrX"],Int64[3e5,2e5,1e4])
genome(chrinfo)
chr_names(chrinfo)
chr_lengths(chrinfo)
chr_ends(chrinfo)
chr_offsets(chrinfo)
chrinfo[2] # 5e5
```
"""
immutable GenomeInfo{T1<:Integer}
    name::String
    chr_ends::AxisArray{T1,1,Vector{Int64},Tuple{AxisArrays.Axis{:chromosome,Array{String,1}}}}
    function GenomeInfo{T1}(name::String, chromosomes::Vector{String}, lengths::Vector{T1})
        length(chromosomes) != length(lengths) && throw(ArgumentError("'chromosomes' and 'lengths' must be the same length."))
        ends = AxisArray(cumsum(lengths), Axis{:chromosome}(chromosomes))
        new(name, ends)
    end
end

function GenomeInfo{T1<:Integer}( name::String, chromosomes::Vector{String}, lengths::Vector{T1})
    GenomeInfo{T1}(name, chromosomes, lengths)
end

function Base.show(io::IO, x::GenomeInfo)
    t = typeof(x)::DataType
    show(io, t)
    write(io,"\nGenome: ", x.name)
    write(io,"\nChromosome Lengths: \n")
    for (name, len) in zip(chr_names(x), chr_lengths(x))
       write(io,"\n $(name) : $(len)")
    end
end

genome(x::GenomeInfo) =  x.name
chr_names(x::GenomeInfo) = axisvalues(x.chr_ends)[1]
chr_ends(x::GenomeInfo) = copy(x.chr_ends)

function chr_offsets(x::GenomeInfo)
    ends = chr_ends(x)
    n = length(ends)
    ends[2:n] = ends[1:(n-1)]
    ends[1] = 0
    ends
end

function chr_lengths(x::GenomeInfo)
    lens = chr_ends(x)
    lens[2:end] = diff(lens)
    lens
end

==(x::GenomeInfo, y::GenomeInfo) = isequal(x.name,y.name) && isequal(x.chr_ends,y.chr_ends)

Base.getindex(x::GenomeInfo, i::Union{String, Integer}) = x.chr_ends[i]

@doc (@doc GenomeInfo) chr_names, chr_lengths, chr_ends, chr_offsets, genome
