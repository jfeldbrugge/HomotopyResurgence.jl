module HomotopyMakieExt

# Makie methods of the plotting functions declared in Homotopy (src/Plots.jl). The extension loads
# with Makie, i.e. as soon as the user loads a backend such as CairoMakie or GLMakie.

using Homotopy, Makie
import Homotopy: plotHomotopy

function plotHomotopy(pathsBorel, K, plotRange = 1, index = length(pathsBorel[1]))
    fig = Figure()
    ax = Axis(fig[1,1], xlabel = "Re", ylabel = "Im", limits=(-plotRange, plotRange, -plotRange, plotRange), aspect=1)
    for (i, path) in pairs(pathsBorel)
        lines!(ax, real.(path), imag.(path), label = string(i))
    end

    points = [path[index] for path in pathsBorel]
    inside(z) = isfinite(z) && sqrt(2) * abs(z) < plotRange

    for i in axes(K, 1), j in axes(K, 2)
        if K[i, j] != 0 && inside(points[i]) && inside(points[j])
            lines!(ax, [Point2f(real(points[i]), imag(points[i])), Point2f(real(points[j]), imag(points[j]))],
                   color = :red, linewidth = 2)
        end
    end

    for i in eachindex(points)
        if inside(points[i])
            scatter!(ax, real(points[i]), imag(points[i]), color = :black)
            text!(ax, string(i), position = Point2f(real(points[i]), imag(points[i])), align = (:right, :bottom), offset = (8, 8))
        end
    end

    axislegend(ax)
    fig, ax
end

end
