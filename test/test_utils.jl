module TestUtils

using GenomicVectors
using Base.Test

@testset begin

    # genopos
    chrinfo = GenomeInfo("hg19",["chr7","chr10","chrM"],Int64[3e4,2e4,1e3])
    g = genopos([1, 405, 123, 7, 12], ["chr10","chr7","chrM","chr10","chr7"], chrinfo)
    @test g == [30001, 405, 50123, 30007, 12]
    @test_throws ArgumentError genopos([1,2,3], ["chr7","chrM"], chrinfo)

    # chrpos
    chrinfo = GenomeInfo("hg19",["chr7","chr10","chrM"],Int64[3e4,2e4,1e3])
    g == [30001, 405, 50123, 30007, 12]
    @test chrpos(g,chrinfo) == [1, 405, 123, 7, 12]

    # chromosomes
    chrinfo = GenomeInfo("hg19",["chr7","chr10","chrM"],Int64[3e4,2e4,1e3])
    g == [30001, 405, 50123, 30007, 12]
    @test chromosomes(g,chrinfo) == ["chr10","chr7","chrM","chr10","chr7"]

    # chr_and_pos
end #testset

end #module
