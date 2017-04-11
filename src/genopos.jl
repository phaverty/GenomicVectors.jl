using ProfileView
using GenomicRanges
using JLD

#function genopos(pos::Vector{Int64}, chromosomes::Vector{String}, info::GenomeInfo{Int64})
function genopos(pos, chromosomes, info)
    offsets = chr_offsets(info)
    lengths = chr_lengths(info)
    gpos = similar(offsets, length(pos))
    chr_hash = info.chr_names
    for i in eachindex(gpos)
        chr = chromosomes[i]
        x = pos[i]
        chr_ind = chr_hash[ chr ]
        if 1 <= x <= lengths[chr_ind]
            gpos[i] = x + offsets[chr_ind]
        else
            error("Position $x is outside the bounds of chromosome $chr (length $(lengths[chr_ind])).")
        end
    end
    gpos
end
