module TestBam

using GenomicVectors
using Base.Test
using BioAlignments

#@testset begin

    bam_file = "GSE25840_GSM424320_GM06985_gencode_spliced.head.bam"
    #bam_path = joinpath(Pkg.dir("BioFmtSpecimens"),"BAM", bam_file)
    # Use a copy of file from BioFmtSpecimens until that package gets registered
    bam_path = joinpath(Pkg.dir("GenomicVectors"),"BAM", bam_file)
    reader = open(BAM.Reader, bam_path)

    info = GenomeInfo("hg19", reader)
    @test isa(info, GenomeInfo)
    gr = GenomicRanges("hg19", reader)
    @test isa(gr, GenomicRanges)

    close(reader)

#end # testset

end # module
