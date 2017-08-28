module TestGenomeInfo

using GenomicVectors
using Base.Test
using GenomicFeatures

@testset begin

    ## GenomeInfo
    chrs = ["chr1","chr2","chrX"]
    x = GenomeInfo("hg19",chrs,Int64[3e5,2e5,1e4])
    @test genome(x) == "hg19"
    @test chr_lengths(x) == Int64[3e5,2e5,1e4]
    @test chr_offsets(x) == Int64[0,3e5,5e5]
    @test chr_ends(x) == cumsum(Int64[3e5,2e5,1e4])
    @test chr_names(x) == chrs
    io = IOBuffer()
    @test typeof(show(io,x)) == Void # At least test that show does not give error
    @test x["chr2"] == 500000
    @test x[2] == 500000

    ## same genome
    chrinfo = GenomeInfo("hg19",["chr1","chr2","chrX"],Int64[3e5,2e5,1e4])
    chrs = ["chr1","chr2","chr2","chrX"]
    s = [100, 200, 300, 400]
    e = [120, 240, 350, 455]
    gr = GenomicRanges(chrs,s,e,chrinfo)

    @test same_genome(gr,Interval("hg19",1,3,'.')) == true
    @test same_genome(Interval("hg19",1,3,'.'),gr) == true
    @test same_genome(Interval("hg18",1,3,'.'),gr) == false
    @test same_genome(Interval("hg19",1,310000,'.'),gr) == false

end # testset

end # module
