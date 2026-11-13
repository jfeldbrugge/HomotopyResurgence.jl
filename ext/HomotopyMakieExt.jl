module HomotopyMakieExt

# Makie methods of the plotting functions declared in Homotopy (src/Plots.jl). The extension loads
# with Makie, i.e. as soon as the user loads a backend such as CairoMakie or GLMakie.

using Homotopy, Makie
import Homotopy: plotHomotopy

function plotHomotopy(pathsBorel, K, plotRange = 1, index = length(pathsBorel[1]); labels = true,
                      legend = length(pathsBorel) ≤ 10)
    fig = Figure()
    ax = Axis(fig[1,1], xlabel = "Re", ylabel = "Im", limits=(-plotRange, plotRange, -plotRange, plotRange), aspect=1)
    for (i, path) in pairs(pathsBorel)
        lines!(ax, real.(path), imag.(path), label = string(i))
    end

    points = [path[index] for path in pathsBorel]

    # edges of the adjacency graph, with a width proportional to |K|; the axis clips them
    for i in axes(K, 1), j in axes(K, 2)
        if i < j && (K[i, j] != 0 || K[j, i] != 0) && isfinite(points[i]) && isfinite(points[j])
            lines!(ax, [Point2f(real(points[i]), imag(points[i])), Point2f(real(points[j]), imag(points[j]))],
                   color = :red, linewidth = 2 * max(abs(K[i, j]), abs(K[j, i])))
        end
    end

    finite = filter(i -> isfinite(points[i]), eachindex(points))
    scatter!(ax, real.(points[finite]), imag.(points[finite]), color = :black)
    if labels
        for i in finite
            abs(real(points[i])) < plotRange && abs(imag(points[i])) < plotRange || continue
            text!(ax, string(i), position = Point2f(real(points[i]), imag(points[i])), align = (:right, :bottom), offset = (8, 8))
        end
    end

    legend && axislegend(ax)
    fig, ax
end

end
