### GRanges Type
type GRanges <: AbstractGRanges
    root::Node
    n::Int
    info::GenomeInfo
    function GRanges(range::IntervalTree, info::GenomeInfo)

    end
end

function GRanges{T1 <: String}(ranges::IntervalTree, chromosomes::Array{T1}, info::GenomeInfo)

end

function GRanges{T1 <: String}(ranges::Array{Interval}, chromosomes::Array{T1}, info::GenomeInfo)

end

function GRanges{T1 <: Integer, T2 <: String}(starts::Array{T1}, ends::Array{T1}, chromosomes::Array{T2}, info::GenomeInfo)

end

function Base.show(io::IO, it::GRanges)
    t = typeof(it)::DataType
    show(io, t)
    n = length(it)
    write(io,"Length: ", n, " intervals")
    if length(it) > 6
        # Hacky random access ...
        for (i, x) in enumerate(it)
            i < 4 && print(io,x)
        end
        write(io,"\n\u22EE") # Vertical ellipsis
        for (i, x) in enumerate(it)
            i > (n-3) && print(io,x)
        end
    else
        for x in it
            write(io,x)
        end
    end
    show(io,x.info)
end
