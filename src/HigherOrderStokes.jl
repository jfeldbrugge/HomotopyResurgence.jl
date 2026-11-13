# Stokes constants follow the convention of IntegralResurgence.jl: for the integral
# (k/2πi)^{d/2} ∫ e^{ik f(x)} dx, the Borel singularity of saddle x_j is t_j = -i f(x_j), and
# K[l, j] = K_{l←j} is the Stokes constant with which the Laplace contour of saddle j crosses t_l.
# The asymptotic series of saddle j is normalised by T0_j = 1/ρ_j with ρ_j² = det Hess f(x_j).
# Stokes constants are bilinear in the normalisations: changing ρ_j → s_j ρ_j (s_j = ±1)
# changes K[l, j] → s_l s_j K[l, j].

"""
    referenceStokes1D(m)

Stokes constants of the reference exponent g(x) = x^m/m - x for the saddles e_n = e^{2πin/(m-1)},
n = 0, …, m-2, with the continuous normalisation ρ_n = √(m-1) e^{-iπn/(m-1)} (ρ_n² = g''(e_n)):
K[l, j] = sign(l - j). Checked against the Gauss-Manin Stokes constants of IntegralResurgence.jl
for m = 3, …, 7. The constants are unchanged when g is multiplied by a phase c, with ρ_n → √c ρ_n.
"""
referenceStokes1D(m) = [sign(l - j) for l in 1:m-1, j in 1:m-1]

"""
    initialSaddles(M)

Saddle points of the factorizable reference system g(x) = Σᵢ cᵢ (xᵢ^mᵢ/mᵢ - xᵢ) (the saddles do
not depend on the phases cᵢ), their signed Stokes constants K and the normalisations ρ of the
series (up to a factor common to all saddles, which does not affect K). Two saddles are adjacent
iff they differ in exactly one coordinate i, with the Stokes constant of the one-dimensional
reference system in that coordinate (`referenceStokes1D`).
"""
function initialSaddles(M)
    labels = vec([collect(t) for t in Iterators.product((0:m - 2 for m in M)...)])
    saddles = [[cis(2π * n / (m - 1)) for (n, m) in zip(a, M)] for a in labels]
    K1 = [referenceStokes1D(m) for m in M]

    K = zeros(Int, length(labels), length(labels))
    for (j, a) in pairs(labels), (l, b) in pairs(labels)
        differ = findall(a .!= b)
        if length(differ) == 1
            i = differ[1]
            K[l, j] = K1[i][b[i] + 1, a[i] + 1]
        end
    end

    ρ = [prod(sqrt(m - 1) * cis(-π * n / (m - 1)) for (n, m) in zip(a, M)) for a in labels]
    return saddles, K, ρ
end

"""
    referenceSystem(M; c = goldenPhases(length(M)))

The reference exponent g(x) = Σᵢ cᵢ (xᵢ^mᵢ/mᵢ - xᵢ). Distinct generic phases cᵢ are essential when
two mᵢ agree: with cᵢ = 1 different saddles share their critical value and many triples are
collinear in the Borel plane, so the reference Stokes constants would be ill defined.
"""
referenceSystem(M; c = goldenPhases(length(M))) = x -> sum(c[i] * (x[i]^M[i] / M[i] - x[i]) for i in eachindex(M))

"Phases cᵢ = e^{iϕ(i-1)} with the golden angle ϕ."
goldenPhases(d) = [cis(π * (3 - sqrt(5)) * (i - 1)) for i in 1:d]

# Indices i where the sign of x changes between x[i] and the next nonzero entry.
# Exact zeros are skipped, so a zero on a grid point gives one crossing (not two),
# and a touch without a sign change gives none.
function zero_crossings(x)
    crossings = Int[]
    last = 0
    for i in eachindex(x)
        (iszero(x[i]) || isnan(x[i])) && continue
        if last != 0 && sign(x[last]) != sign(x[i])
            push!(crossings, i - 1)
        end
        last = i
    end
    return crossings
end

function findColinear(z1, z2, z3)
    δλ = 1 / (length(z1) - 1)
    volume = imag((z2 .- z1) .* conj.(z3 - z1))
    crossings = zero_crossings(volume)
    return crossings, δλ .* crossings
end

"Signed area of t_k relative to the directed segment t_j → t_l: > 0 when t_k lies on its left."
side(tj, tk, tl) = imag(conj(tl - tj) * (tk - tj))

"""
Higher-order Stokes events: instants (between the grid points i and i+1) at which the Borel
singularity t_k crosses the segment between t_j and t_l. Returns named tuples (τ, j, k, l) with τ
the interpolated crossing time, sorted by τ.
"""
function hosEvents(pathsBorel)
    n = length(pathsBorel)
    triples = collect(combinations(1:n, 3))
    found = [Vector{NamedTuple{(:τ, :i, :j, :k, :l), Tuple{Float64, Int, Int, Int, Int}}}() for _ in 1:Threads.maxthreadid()]
    nonzero(x) = !iszero(x) && !isnan(x)
    Threads.@threads for tr in triples
        a, b, c = pathsBorel[tr[1]], pathsBorel[tr[2]], pathsBorel[tr[3]]
        any(isnan, (a[1], b[1], c[1])) && continue
        area = imag.((b .- a) .* conj.(c .- a))
        for cr in zero_crossings(area)
            # the sign changes between the nonzero areas at i1 and i2 (exact zeros in between are skipped)
            i2 = cr + 1
            i1 = findprev(nonzero, area, cr)
            τ = i1 + area[i1] / (area[i1] - area[i2]) * (i2 - i1)
            w = τ - floor(τ); i0 = floor(Int, τ)
            z = [(1 - w) * p[i0] + w * p[min(i0 + 1, length(p))] for p in (a, b, c)]
            # the middle point lies inside the segment of the other two
            for (m, e1, e2) in ((2, 1, 3), (1, 2, 3), (3, 1, 2))
                s = real(conj(z[e2] - z[e1]) * (z[m] - z[e1])) / abs2(z[e2] - z[e1])
                if 0 < s < 1
                    push!(found[Threads.threadid()], (τ = τ, i = i1, j = tr[e1], k = tr[m], l = tr[e2]))
                    break
                end
            end
        end
    end
    sort(reduce(vcat, found), by = e -> e.τ)
end

"""
    higherOrderStokes(K, pathsBorel)

Continue the Stokes constants K (K[l, j] = K_{l←j}, normalisations continuous along the paths)
through the higher-order Stokes phenomena of the homotopy. When t_k crosses the segment from t_j
to t_l from its left to its right,
    K_{l←j} ↦ K_{l←j} + K_{k←j} K_{l←k},
and with the opposite sign when it crosses from right to left (the Laplace contour of j passes
t_k on the other side; the difference is the Stokes jump of j at t_k, which then jumps at t_l).
The same event updates K_{j←l} with the roles of j and l exchanged.
"""
function higherOrderStokes(K, pathsBorel)
    K = copy(K)
    for e in hosEvents(pathsBorel)
        i = e.i
        for (j, l) in ((e.j, e.l), (e.l, e.j))
            p = pathsBorel
            before = side(p[j][i], p[e.k][i], p[l][i])
            σ = before > 0 ? 1 : -1                       # left → right: +, right → left: -
            K[l, j] += σ * K[e.k, j] * K[l, e.k]
        end
    end
    return K
end

"""
    stokesHomotopy(f, M; c = goldenPhases(length(M)), γ = DEFAULT_γ, N = 10000, escape = 1e6)

Stokes constants of the saddle points of the exponent f (integral ∫ e^{ik f(x)} dx) by homotopy
continuation from the reference system g = Σᵢ cᵢ (xᵢ^mᵢ/mᵢ - xᵢ), with mᵢ - 1 ≥ deg ∂f/∂xᵢ. The
saddles of g move along h(x, λ) = λ f(x) + γ(1 - λ) g(x) into the saddles of f or to infinity,
and the Stokes constants change at every higher-order Stokes phenomenon on the way.

Returns a named tuple with the saddle points `x` of f that are reached, their Borel singularities
`Σ = -i f(x)`, the Stokes constants `K` (K[l, j] = K_{l←j}, as `stokes_constants` of
IntegralResurgence.jl) for the normalisations `T0` of the series (see `seriesNormalisation`), and
the tracked `paths`, `pathsBorel` and continuous normalisations `ρs` of all reference saddles, with
`finite` the indices of the paths that reach f.
"""
function stokesHomotopy(f, M; c = goldenPhases(length(M)), γ = DEFAULT_γ, N = 10000, escape = 1e6,
                        normalisation = :morse)
    g = referenceSystem(M; c = c)
    saddles, K0, ρ0 = initialSaddles(M)
    # ρ² = det Hess(γ g) = γ^d Πᵢ cᵢ g''; the common factor √(γ^d Πᵢ cᵢ) does not affect K
    ρ0 = sqrt(γ)^length(M) * prod(sqrt.(c)) .* ρ0
    paths, pathsBorel, ρs = followHomotopy(saddles, [f, g], N, γ; ρ0 = ρ0, escape = escape)
    K = higherOrderStokes(K0, pathsBorel)

    finite = [i for i in eachindex(paths) if all(isfinite, paths[i][end])]
    x = [refineSaddle(f, paths[i][end], 1e-13)[1] for i in finite]
    Σ = [-im * f(p) for p in x]
    # from the continuous branch ρ along the path to the normalisation of the series of f
    ρf = [seriesNormalisation(cHessian(f, p), p; normalisation) for p in x]
    s = [round(Int, real(ρf[a] / ρs[i][end])) for (a, i) in pairs(finite)]
    Kf = [s[a] * s[b] * K[finite[a], finite[b]] for a in eachindex(finite), b in eachindex(finite)]

    return (; x, Σ, K = Kf, T0 = 1 ./ ρf, paths, pathsBorel, ρs, finite, Kfull = K)
end

"""
    seriesNormalisation(H, x; normalisation = :morse, atol = 1e-8)

The branch ρ of √det H that normalises the asymptotic series of the saddle point x with Hessian
H, T0 = 1/ρ. Stokes constants are only defined once these branches are fixed.

- `:morse` (default, as `transseries` of IntegralResurgence.jl): the principal branch for a
  complex saddle point and, for a real saddle point with Morse index ν, ρ = i^ν √|det H|, the
  stationary-phase normalisation.
- `:principal`: the principal branch. When det H lies on the negative real axis (up to `atol`),
  the branch is fixed at ρ = +i √|det H|; otherwise the sign would depend on the rounding error
  in the imaginary part of det H.
"""
function seriesNormalisation(H, x; normalisation = :morse, atol = 1e-8)
    d = det(H)
    scale = max(1, maximum(abs, H))
    isreal_point = maximum(abs ∘ imag, x) < atol * max(1, maximum(abs, x)) &&
                   maximum(abs ∘ imag, H) < atol * scale
    if normalisation == :morse && isreal_point
        ν = count(<(0), eigvals(Symmetric(real.(H))))
        return im^ν * sqrt(abs(d))
    elseif normalisation in (:morse, :principal)
        abs(imag(d)) < atol * abs(d) && real(d) < 0 && return im * sqrt(abs(d))
        return sqrt(complex(d))
    end
    error("seriesNormalisation: normalisation must be :morse or :principal")
end
