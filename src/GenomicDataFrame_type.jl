const ColumnIndex = Union{Symbol,Integer}

"""
# GenomicDataFrame

A DataFrame-like class with a GenomicVector as an index.

# Examples
```julia
    genomeinfo = GenomeInfo("hg19",["chr1","chr2","chrX"],Int64[3e5,2e5,1e4])
    chrs = ["chr2","chr1","chr2","chrX"]
    pos = Int64[3e4,4.2e3,1.9e5,1e4]
    gp = GenomicPositions(pos,chrs,genomeinfo)
    dt = DataFrame(a=1:4,b=5:8)
    gt = GenomicDataFrame(gp,dt)
    gt[1:2,1:2]
```
"""
struct GenomicDataFrame{T1 <: AbstractGenomicVector, T2 <: AbstractDataFrame} <: AbstractDataFrame
    rowindex::T1
    table::T2
function GenomicDataFrame{T1,T2}(rowindex,table) where {T1 <: AbstractGenomicVector,T2 <: AbstractDataFrame}
    nrow = size(table,1)
        if nrow != 0 && length(rowindex) != nrow
            throw(ArgumentError("GenomicDataFrame requires that `length(rowindex) == size(table,1)`"))
        end
        new(rowindex,table)
    end
end
GenomicDataFrame(gv,dt) = GenomicDataFrame{typeof(gv),typeof(dt)}(gv,dt)
rowindex(gt::GenomicDataFrame) = copy(gt.rowindex)
table(gt::GenomicDataFrame) = copy(gt.table)
_rowindex(gt::GenomicDataFrame) = gt.rowindex
_table(gt::GenomicDataFrame) = gt.table
for op in [:size, :nrow, :ncol, :index, :columns]
    @eval (DataFrames.$op)(x::GenomicDataFrame) = ($op)(_table(x))
end
Base.names(x::GenomicDataFrame) = names(_table(x))

function Base.show(io::IO, x::GenomicDataFrame)
    t = typeof(x)
    show(io, t)
    println("\n\nRow Index:")
    show(io,_rowindex(x))
    println("\n\nTable:")
    show(io,_table(x))
end

Base.getindex(gt::GenomicDataFrame,i,j) = GenomicDataFrame(rowindex(gt)[i], table(gt)[i,j])
Base.getindex(gt::GenomicDataFrame,j) = GenomicDataFrame(rowindex(gt), table(gt)[j])
Base.getindex(gt::GenomicDataFrame,j::ColumnIndex) = table(gt)[j]
Base.setindex!(gt::GenomicDataFrame,value,j) = setindex!(_table(gt),value,j)
Base.setindex!(gt::GenomicDataFrame,value,i,j) = setindex!(_table(gt),value,i,j)

## Getters that delegate to the AbstractDataFrame row index
for op in [:chr_info, :_strands, :_genostarts, :_genoends, :starts, :ends, :widths, :chromosomes, :genostarts, :genoends, :strands, :each, :chrpos, :genopos, :chrindex]
    @eval ($op)(x::GenomicDataFrame) = ($op)(_rowindex(x))
end
## Searches that delegate to the genome info
for op in [:indexin, :findin, :in]
    @eval (Base.$op)(x::GenomicDataFrame,y::AbstractGenomicVector,exact::Bool=true) = ($op)(_rowindex(x),y,exact)
    @eval (Base.$op)(x::AbstractGenomicVector,y::GenomicDataFrame,exact::Bool=true) = ($op)(x,_rowindex(y),exact)
    @eval (Base.$op)(x::GenomicDataFrame,y::GenomicDataFrame,exact::Bool=true) = ($op)(_rowindex(x),_rowindex(y),exact)
end

## Table ops that delegate to the table
#for op in [:join]
#
#end

function Base.vcat(x::GenomicDataFrame,y::GenomicDataFrame)
    same_genome(x,y) || throw(ArgumentError("x and y must be from the same genome"))
    GenomicDataFrame( vcat(_rowindex(x),_rowindex(y)), vcat(_table(x),_table(y)) )
end
