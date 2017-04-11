### Extra methods for IntervalTree.jl

## Getter for IntervalValue
value{K, V}(i::IntervalValue{K, V}) = i.value

## Display for Interval, IntervalValue, IntervalTree
Base.print(io::IO, x::Interval) = print(io, "\n($(first(x)),$(last(x)))")
function Base.show(io::IO, x::Interval)
    t = typeof(x)::DataType
    show(io, t)
    print(x)
end

Base.print(io::IO, x::IntervalValue) = print(io, "\n($(first(x)),$(last(x))) => $(value(x))")
function Base.show(io::IO, x::IntervalValue)
    t = typeof(x)::DataType
    show(io, t)
    print(x)
end

function Base.show(io::IO, it::IntervalTree)
    t = typeof(it)::DataType
    show(io, t)
    n = length(it)
    println("Length: ", n, " intervals")
    if length(it) > 6
        # Hacky random access ...
        for (i, x) in enumerate(it)
            i < 4 && print(x)
        end
        print("\n\u22EE") # Vertical ellipsis
        for (i, x) in enumerate(it)
            i > (n-3) && print(x)
        end
    else
        for x in it
            print(x)
        end
    end
end
