#########################
### Utility Functions ###
#########################

"""
Given chromosome and chromosome position information and a description of
the chromosomes (a GenoPos object), calculate the corresponding positions
in the linear genome.
"""
function genopos(pos, chromosomes, chrinfo::GenomeInfo)
    offsets = chr_offsets(chrinfo)
    lengths = chr_lengths(chrinfo)
    gpos = similar(offsets, length(pos))
    prev_chr = chromosomes[1]
    len = lengths[prev_chr]
    @inbounds for (i,x,chr) in zip(1:length(pos), pos, chromosomes)
        if chr != prev_chr
            prev_chr = chr
            len = lengths[prev_chr]
        end
        if 1 <= x <= len
            gpos[i] = x + offsets[prev_chr]
        else
            error("Position $x is outside the bounds of chromosome $chr (length $(lengths[prev_chr])).")
        end
    end
    gpos
end

"""
Given positions in the linear genome, calculate the position on the relevant chromosome.
"""
function chrpos(pos, chrinfo::GenomeInfo)
    ends = chr_ends(chrinfo)
    offsets = chr_offsets(chrinfo)
    nchr = length(ends)
    res = similar(pos,length(pos))
    r = 1
    i = 1
    @inbounds for g in pos
        if g > ends[r] || g <= offsets[r]
            r = searchsortedfirst(ends, g, 1, nchr, Base.Forward)
        end
        res[i] = pos[i] - offsets[ r ]
        i = i + 1
    end
    res
end
