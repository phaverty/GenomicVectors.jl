#########################
### Utility Functions ###
#########################

"""
Given chromosome and chromosome position information and a description of
the chromosomes (a GenoPos object), calculate the corresponding positions
in the linear genome.
"""
function genopos(positions, chromosomes, chrinfo::GenomeInfo)
    if length(positions) != length(chromosomes)
        throw(ArgumentError("Arguments positions and chromosomes must have the same length."))
    end
    offsets = chr_offsets(chrinfo)
    lengths = chr_lengths(chrinfo)
    gpos = similar(offsets, length(positions))
    prev_chr = chromosomes[1]
    len = lengths[prev_chr]
    @inbounds for (i,x,chr) in zip(1:length(positions), positions, chromosomes)
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
function chrpos(positions, chrinfo::GenomeInfo)
    ends = chr_ends(chrinfo)
    offsets = chr_offsets(chrinfo)
    nchr = length(ends)
    res = similar(positions,length(positions))
    r = 1
    i = 1
    @inbounds for g in positions
        if g > ends[r] || g <= offsets[r]
            r = searchsortedfirst(ends, g, 1, nchr, Base.Forward)
        end
        res[i] = positions[i] - offsets[ r ]
        i = i + 1
    end
    res
end

"""
Given positions in the linear genome, calculate the position on the relevant chromosome.
"""
function chromosomes(positions, chrinfo::GenomeInfo)
    ends = chr_ends(chrinfo)
    offsets = chr_offsets(chrinfo)
    chrs = chr_names(chrinfo)
    nchr = length(ends)
    res = similar(chr_names(chrinfo),length(positions))
    r = 1
    i = 1
    @inbounds for g in positions
        if g > ends[r] || g <= offsets[r]
            r = searchsortedfirst(ends, g, 1, nchr, Base.Forward)
        end
        res[i] = chrs[r]
        i = i + 1
    end
    res
end

# chr_and_pos (vec of tuples)
