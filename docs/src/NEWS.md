## Version 0.0.1: Initial Public Release
This is the first public release of GenomicVectors.jl, which brings the features of
BioConductor's GenomicRanges to Julia.

## Version 0.0.2: changes to indexing a GenomicRanges
This release changes how indexing a GenomicRanges at a single position works. From the docs:

Indexing a GenomicRanges with an array produces a new GenomicRanges. Getting/setting by a
scalar gives/takes a three-tuple of the (start,end,strand) in genome location units. In the
near future, this will likely switch to using Bio.Intervals.Interval. The each function
produces an iterator of (start,end) two-tuples in genome location units. This is use for many
internal functions, like sorting. This is intentionally similar to RLEVectors.each.

This release also has simplified sorting code, fixes for a few edge cases and 100 percent
test coverage (well, not counting the precompiles).

## Version 0.0.3: Update to scalar indexing, behaves like Vector{Interval}
This version changes the value received or returned by setindex and getindex on a GenomicRanges.
We now use the Interval type from Bio.jl to represent a single range. However, we still maintain
our genome-level coordinates and use the name of the genome (e.g. hg19) in the Interval seqname
slot. Using Interval simplifies some code, but more importantly gives the single range some
meaning and guaranteed behavior. In effect, GenomicRanges now behaves as if it were
Vector{Interval}, but the internal representation is different. Additionally, the metadata slot
of the Interval returned by getindex contains the index of that Interval in the GenomicRanges.
This is useful for recovering the original ordering when a round trip is made converting to
IntervalCollection and back.

The each method on GenomicRanges still returns an iterator of (start,end) tuples, in units of
genomic position, which is useful for sorting and some other internal operations.
