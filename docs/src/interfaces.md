## Interfaces

Much of the functionality of the types in `GenomicVectors.jl` is provided by two
interfaces: the GenomeInfo Interface and the Genomic Location Interface.

### GenomeInfo Interface
This interface provides information regarding the underlying genome of a collection
of locations. Collections from the identical genome may be compared, combined, etc. This
Interface requires that a type implement the following methods:

- chr_info: Returns a `GenomeInfo` object

In return types will be able to provide the GenomeInfo API:

- same_genome: Are two objects of the same genome?
- chr_names: What are the names of genome's chromosomes?
- chr_lengths: What are the lengths of genome's chromosomes?
- chr_ends: What are the endpoint positions of genome's chromosomes?
- chr_offsets: What are the position offsets into each of genome's chromosomes?
- genome: What is the name of the genome?

### Genopos Interface
This interface provides access to the specific locations, relative to a linearized genome
or relative to the chromosome in question. Locations may be queried, moved, size-adusted
and subject to Set operations. For the purposes of comparison of locations the standard
search api is implemented (`findin`, `indexin`, etc.). These all perform exact matching
by default. For each implemented method, an extra argument, `exact` may be set to false
to allow for any type of range overlap to count as a match. (Specific types of overlap,
e.g. complete overlap, overlap of the left bound, etc. may be implemented in the
future.)

This interface requires that the a type implement methods on the following generics:

- _genostarts: the starting nucleotide index for each range/position in the linearized genome
- _genends: the ending nucleotide index for each range/position in the linearized genome
- _strands: the DNA strand for each range/position
- same_genome: tests two

These internal functions are assumed to pass their values by reference. Calling functions
are expected **not** to modify these values.

In turn, the interface provides the following methods:

- starts: the starting nucleotide index for each range/position relative to the chromosome on which they lie
- ends: the ending nucleotide index for each range/position relative to the chromosome on which they lie
- widths: the distance, between the start and end nucleotide of the range, 1 for positions
- chromosomes: the name of the chromosome for each range/position
- genostarts: the starting nucleotide index for each range/position in the linearized genome
- genoends: the ending nucleotide index for each range/position in the linearized genome
- strands: the DNA strand for each range/position, pass by copy
- each: an iterator that returns tuples of genostart and genoend pairs

All of these exported functions should be non-aliasing, or pass-by-copy.
