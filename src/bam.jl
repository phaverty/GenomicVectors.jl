############################################
### Functions for working with BAM files ###
############################################

function GenomeInfo(genome_name, reader::BioAlignments.BAM.Reader)
    GenomeInfo(genome_name,
               reader.refseqnames,
               reader.refseqlens
               )
end

function GenomicRanges(genome_name, reader::BioAlignments.BAM.Reader)
    info = GenomeInfo(genome_name, reader)
    chr = String[]
    left_pos = Int64[]
    right_pos = Int64[]
    strand_pos = Vector{Strand}(0)
    record = BAM.Record()
    ## FIXME: convert to genopos in the loop?
    while !eof(reader)
        read!(reader, record)
        if BAM.ismapped(record)
            push!(chr, BAM.refname(record))
            push!(left_pos, leftposition(record))
            push!(right_pos, rightposition(record))
            push!(strand_pos, strand(record))
        end
    end
    GenomicRanges(chr, left_pos, right_pos, strand_pos, info)
end

function strand(record::BAM.Record)
    ## FIXME: push up to BioAligments
    if BAM.flag(record) & 0x10 == 0
        s = STRAND_POS
    else
        s = STRAND_NEG
    end
    s
end
