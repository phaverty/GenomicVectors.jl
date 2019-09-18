using .RCall

function RCall.sexp(f::GenomeInfo)
    RCall.R"library(GenomicRanges)"
    s = [string(x) for x in chr_names(f)]
    l = collect(chr_lengths(f))
    g = genome(f)
    RCall.R"Seqinfo($s, $l, NA, $g)"
end

function RCall.rcopy(::Type{GenomeInfo}, s::Ptr{RCall.S4Sxp})
    GenomeInfo(
        RCall.rcopy(String, s[:genome][1]),
        RCall.rcopy(Vector{String}, s[:seqnames]),
        RCall.rcopy(Vector{Int64}, s[:seqlengths])
    )
end

function RCall.sexp(f::GenomicRanges)
    RCall.R"library(GenomicRanges)"
    s = starts(f)
    e = ends(f)
    c = [string(x) for x in chromosomes(f)]
    r = convert(Vector{String}, strands(f))
    i = RCall.rcopy(chr_info(f))
#    RCall.R"GRanges($c, IRanges($s, $e), $r, seqinfo = $i)"
#    RCall.R"GRanges($c, IRanges($s, $e), seqinfo = $i)
#    RCall.R"GRanges($c, IRanges($s, $e), $r)
    RCall.R"GRanges($c, IRanges($s, $e))
end

#function RCall.sexp(f::GenomicDataFrame)
#    RCall.R"library(GenomicRanges)"
#    s = starts(f)
#    e = ends(f)
#    c = [string(x) for x in chromosomes(f)]
#    r = convert(Vector{String}, strands(f))
#    i = RCall.rcopy(chr_info(f))
#    m = table(f)
#    RCall.R"GRanges($c, IRanges($s, $e), $r, seqinfo = $i, mcols = DataFrame($m))"
#end

#function RCall.rcopy(::Type{GenomicDataFrame}, s::Ptr{RCall.S4Sxp})
#    GenomicDataFrame(
#        GenomicRanges(
#            RCall.rcopy(Vector{String}, RCall.R"as.character(seqnames($s))"),
#            RCall.rcopy(Vector{Int64},  RCall.R"start($s)"),
#            RCall.rcopy(Vector{Int64},  RCall.R"end($s)"),
#            RCall.rcopy(GenomeInfo, s[:seqinfo])
#        ),
#        RCall.rcopy(DataFrame, RCall.R"as.data.frame(mcols($s))")
#    )
#end

RCall.rcopytype(::Type{RCall.RClass{:Seqinfo}}, s::Ptr{RCall.S4Sxp}) = GenomeInfo
#RCall.rcopytype(::Type{RCall.RClass{:GRanges}}, s::Ptr{RCall.S4Sxp}) = GenomicDataFrame