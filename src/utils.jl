#########################
### Utility Functions ###
#########################

"""
# The GenoPos Interface

Provides access to positional information in the linearized
genome or in chromosome coordinate (e.g. chr4:1000-1020). This interface requires
a type to implement the a method on the non-copying accessors `_genostarts`,
`_genoends` and `_strands` as well as the GenomeInfo Interface.

    starts(x)
Get the starting nucleotide index for each range/position relative to the chromosome on which they lie.

    ends(x)
Get the ending nucleotide index for each range/position relative to the chromosome on which they lie.

    widths(x)
Get the distance, between the start and end nucleotide of the range, An `RleVector` of 1s for single-nucleotide positions.

    chromosomes(x)
    chromosomes(genopos, chrinfo)
Get the name of the chromosome for each range/position.

    genostarts(x)
Get the starting nucleotide index for each range/position in the linearized genome.
    
    genoends(x)
Get the ending nucleotide index for each range/position in the linearized genome.

    strands(x)
Get the DNA strand for each range/position, pass by copy.

    each(x)
Return an iterator that returns tuples of genostart and genoend pairs.

## Utility Functions
Much of this functionality is derived from a few utility functions:

    chrpos(genopos, chrinfo)
Given positions in the linear genome, calculate the position on the relevant chromosome.

    chrindex(genopos, chrinfo)
Given positions in the linear genome, determine the corresponding chromosomes and return the indices of the chromosome in chrinfo.

    genopos(chrpos, chromosomes, chrinfo)
Given chromosome and chromosome position information and a description of
the chromosomes (a GenoPos object), calculate the corresponding positions
in the linear genome.
"""
starts, ends, widths, chromosomes, genostarts, genoends, strands, each, chrpos, genopos, chrindex

function genopos(positions, chromosomes, chrinfo::GenomeInfo)
    if length(positions) != length(chromosomes)
        throw(ArgumentError("Arguments positions and chromosomes must have the same length."))
    end
    offsets = chr_offsets(chrinfo)
    lengths = chr_lengths(chrinfo)
    gpos = similar(offsets, length(positions))
    prev_chr = chromosomes[1]
    len = lengths[prev_chr]
    o = offsets[prev_chr]
    @inbounds for (i,x,chr) in zip(1:length(positions), positions, chromosomes)
        if chr != prev_chr
            prev_chr = chr
            len = lengths[prev_chr]
            o = offsets[prev_chr]
        end
        if 1 <= x <= len
            gpos[i] = x + o
        else
            error("Position $x is outside the bounds of chromosome $chr (length $(lengths[prev_chr])).")
        end
    end
    gpos
end

function chrindex(positions, chrinfo::GenomeInfo)
    ends = chr_ends(chrinfo)
    offsets = chr_offsets(chrinfo)
    res = Vector{Int64}(length(positions))
    i = r = 1
    e = ends[r]
    o = offsets[r]
    @inbounds for g in positions
        if g > e || g <= o
            r = 1
            while g > ends[r]
                r = r + 1
            end
            e = ends[r]
            o = offsets[r]
        end
        res[i] = r
        i = i + 1
    end
    res
end

function chrpos(positions, chrinfo::GenomeInfo)
    ends = chr_ends(chrinfo)
    offsets = chr_offsets(chrinfo)
    res = similar(positions,length(positions))
    i = r = 1
    e = ends[r]
    o = offsets[r]
    @inbounds for g in positions
        if g > e || g <= o
            r = 1
            while g > ends[r]
                r = r + 1
            end
            e = ends[r]
            o = offsets[r]
        end
        res[i] = g - o
        i = i + 1
    end
    res
end

function chromosomes(positions, chrinfo::GenomeInfo)
    ends = chr_ends(chrinfo)
    offsets = chr_offsets(chrinfo)
    chrs = chr_names(chrinfo)
    res = similar(chrs,length(positions))
    i = r = 1
    e = ends[r]
    o = offsets[r]
    @inbounds for g in positions
        if g > ends[r] || g <= offsets[r]
            r = 1
            while g > ends[r]
                r = r + 1
            end
            e = ends[r]
            o = offsets[r]
        end
        res[i] = chrs[r]
        i = i + 1
    end
    res
end

genostarts(x) = copy(_genostarts(x))
genoends(x) = copy(_genoends(x))
strands(x) = copy(_strands(x))
RLEVectors.starts(x) = chrpos(_genostarts(x),chr_info(x))
RLEVectors.ends(x) = chrpos(_genoends(x),chr_info(x))
RLEVectors.widths(x) = (_genoends(x) - _genostarts(x)) + 1
RLEVectors.each(x) = zip(_genostarts(x),_genoends(x))
chromosomes(x) = chromosomes(_genostarts(x),chr_info(x))

# Other candidates for GenoPos Interface or AbstractGenomicVector include Vector{Interval}, issorted, sortperm, show, convert(DataTable,x)
