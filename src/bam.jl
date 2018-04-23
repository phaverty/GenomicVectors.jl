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
    for record in reader
        push!(chr, BAM.refname(record))
        push!(left_pos, leftposition(record))
        push!(right_pos, rightposition(record))
    end
    GenomicRanges(chr, left_pos, right_pos, info)
end
