typealias ColumnIndex Union{Symbol,Integer}

"""A DataTable-like class with a GenomicVector as an index

### Constructors
```julia
    genomeinfo = GenomeInfo("hg19",["chr1","chr2","chrX"],Int64[3e5,2e5,1e4])
    chrs = ["chr2","chr1","chr2","chrX"]
    pos = Int64[3e4,4.2e3,1.9e5,1e4]
    gp = GenomicPositions(pos,chrs,genomeinfo)
    dt = DataTable(a=1:4,b=5:8)
    gt = GenomicTable(gp,dt)
    gt[1:2,1:2]
```
"""
type GenomicTable{T1 <: AbstractGenomicVector, T2 <: AbstractDataTable} <: AbstractDataTable
    rowindex::T1
    table::T2
    function GenomicTable(rowindex,table)
        if length(rowindex) != nrow(table)
            throw(ArgumentError("GenomicTable requires that `length(rowindex) == nrow(table)`"))
        end
        new(rowindex,table)
    end
end
GenomicTable(gv,dt) = GenomicTable{typeof(gv),typeof(dt)}(gv,dt)
rowindex(gt::GenomicTable) = copy(gt.rowindex)
table(gt::GenomicTable) = copy(gt.table)
_rowindex(gt::GenomicTable) = gt.rowindex
_table(gt::GenomicTable) = gt.table

function Base.show(io::IO, x::GenomicTable)
    t = typeof(x)
    show(io, t)
    println("\n\nRow Index:")
    show(io,_rowindex(x))
    println("\n\nTable:")
    show(io,_table(x)
end

Base.getindex(gt::GenomicTable,i,j) = GenomicTable(rowindex(gt)[i], table(gt)[i,j])
Base.getindex(gt::GenomicTable,j) = GenomicTable(rowindex(gt), table(gt)[j])
Base.getindex(gt::GenomicTable,j::ColumnIndex) = table(gt)[j]
