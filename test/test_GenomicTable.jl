module TestGenomicPositions

using GenomicVectors
using Base.Test
using DataTables

@testset begin

    genomeinfo = GenomeInfo("hg19",["chr1","chr2","chrX"],Int64[3e5,2e5,1e4])
    chrs = ["chr2","chr1","chr2","chrX"]
    pos = Int64[3e4,4.2e3,1.9e5,1e4]
    gp = GenomicPositions(pos,chrs,genomeinfo)
    dt = DataTable(a=1:4,b=5:8)
    gt = GenomicTable(gp,dt)
    @test isa(gt,GenomicTable)
    @test size(gt) == (4,2)
    gt = gt[1:2,1:2]
    @test isa(gt,GenomicTable)
    @test size(gt) == (2,2)
    @test table(gt) == GenomicVectors._table(gt)
    @test rowindex(gt) == GenomicVectors._rowindex(gt)
    @test table(gt) == DataTable(a=1:2,b=5:6)
    @test rowindex(gt) == GenomicPositions(pos[1:2],chrs[1:2],genomeinfo)
    gt = GenomicTable(gp,dt)
    @test convert(Vector,gt[:b]) == convert(Vector,gt[2])
    gt[2] = [4,3,2,1]
    @test convert(Vector,gt[2]) == [4,3,2,1]
    gt[2:3,[2]] = [7,6]
    @test convert(Vector,gt[2]) == [4,7,6,1]
    
end # testset

end # module
