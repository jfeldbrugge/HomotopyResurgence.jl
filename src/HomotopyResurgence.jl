module HomotopyResurgence
    using LinearAlgebra, ForwardDiff, Combinatorics, Base.Threads

    include("HomotopyContinuation.jl")
    include("HigherOrderStokes.jl")
    include("Plots.jl")

    export homotopy, followHomotopy, findColinear, higherOrderStokes, plotHomotopy, initialSaddles
    export stokesHomotopy, referenceSystem, referenceStokes1D, goldenPhases, seriesNormalisation

    # Calling a plotting function without a Makie backend gives a MethodError; point to the extension.
    function __init__()
        Base.Experimental.register_error_hint(MethodError) do io, exc, argtypes, kwargs
            if exc.f isa Function && nameof(exc.f) in MAKIE_FUNCTIONS && isnothing(Base.get_extension(HomotopyResurgence, :HomotopyResurgenceMakieExt))
                print(io, "\n\n", nameof(exc.f), " needs a Makie backend. Load one first, e.g. `using CairoMakie`.")
            end
        end
    end
end
