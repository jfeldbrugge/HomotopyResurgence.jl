# Plotting. The methods live in the package extension ext/HomotopyResurgenceMakieExt.jl, which loads with
# Makie: run `using CairoMakie` (or GLMakie) to enable them. Without a Makie backend the functions
# exist but have no methods, and calling them prints a hint (see `__init__`).

"""
    plotHomotopy(pathsBorel, K, plotRange = 1, index = length(pathsBorel[1]))

Plot the tracks of the Borel singularities in the homotopy together with the adjacency graph
(an edge for every K[i, j] ≠ 0) at the grid point `index`. Returns `(fig, ax)`.
Requires a Makie backend: load `CairoMakie` (or `GLMakie`) first.
"""
function plotHomotopy end

const MAKIE_FUNCTIONS = (:plotHomotopy,)
