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

    genopos(chrpos, chromosomes, chrinfo)
Given chromosome and chromosome position information and a description of
the chromosomes (a GenoPos object), calculate the corresponding positions
in the linear genome.
"""
starts, ends, widths, chromosomes, genostarts, genoends, strands, each, chrpos, genopos

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
        r = min(r,nchr)
        res[i] = positions[i] - offsets[ r ]
        i = i + 1
    end
    res
end

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
        r = min(r,nchr)
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
