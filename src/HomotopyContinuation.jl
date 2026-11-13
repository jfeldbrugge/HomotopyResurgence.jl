function cGradient(f, p)
    Z = zeros(length(p))
    R = ForwardDiff.gradient(pp -> real(f(pp + p)), Z)
    I = ForwardDiff.gradient(pp -> imag(f(pp + p)), Z)
    return R + im * I
end

function cHessian(f, p)
    Z = zeros(length(p))
    R = ForwardDiff.hessian(pp -> real(f(pp + p)), Z)
    I = ForwardDiff.hessian(pp -> imag(f(pp + p)), Z)
    return R + im * I
end

function refineSaddle(f, x0, tol = 1e-5, max_iterations = 100)
    x = copy(x0)
    ∇f = cGradient(f, x)
    res_norm = norm(∇f)

    for k in 1:max_iterations
        if res_norm ≤tol
            return x, res_norm
        else
            J = cHessian(f, x)
            Δ = -(J \ ∇f)
            x = x + Δ

            ∇f = cGradient(f, x)
            res_norm = norm(∇f)
        end
    end
    return x, res_norm
end

# Generic complex phase for the γ-trick; a real γ (such as 1) lets paths collide.
const DEFAULT_γ = cis(1.2345)

"""
One predictor-corrector step of the saddle point p of h(⋅, λ) from λa to λb. The step is halved
until Newton's method converges quickly with a small correction, so that the path cannot jump to
a neighbouring saddle point. Returns `nothing` when the step cannot be completed.
"""
function trackStep(f, g, γ, p, λa, λb; tol = 1e-10, depth = 0, max_depth = 30)
    Hλ(λ) = λ .* cHessian(f, p) .+ γ * (1 - λ) .* cHessian(g, p)
    Δ = -(Hλ(λa) \ (cGradient(f, p) .- γ .* cGradient(g, p)))
    q = p .+ (λb - λa) .* Δ
    h(x) = λb * f(x) + γ * (1 - λb) * g(x)
    scale = 1 + norm(p)
    converged = false
    for _ in 1:8
        G = cGradient(h, q)
        all(isfinite, G) || break
        δ = -(cHessian(h, q) \ G)
        q = q .+ δ
        if norm(δ) ≤ tol * (1 + norm(q))
            converged = true
            break
        end
    end
    if converged && norm(q .- p) ≤ 0.1 * scale
        return q
    end
    depth ≥ max_depth && return nothing
    λm = (λa + λb) / 2
    pm = trackStep(f, g, γ, p, λa, λm; tol, depth = depth + 1, max_depth)
    pm === nothing && return nothing
    return trackStep(f, g, γ, pm, λm, λb; tol, depth = depth + 1, max_depth)
end

"""
    homotopy(p, Fs, N, γ = DEFAULT_γ; escape = 1e6)

Track the saddle point p of γ g (with Fs = [f, g]) along h(x, λ) = λ f(x) + γ (1 - λ) g(x) on the
grid λ = 0, 1/N, ..., 1, with adaptive substeps between grid points. A path that leaves the ball
|x| < escape, or cannot be continued, goes to infinity: its remaining points are NaN.
"""
function homotopy(p, Fs, N, γ = DEFAULT_γ; escape = 1e6)
    f, g = Fs[1], Fs[2]
    λs = collect(range(0, 1, length = N + 1))
    path = Vector{Vector{ComplexF64}}(undef, N + 1)
    path[1] = ComplexF64.(p)
    nanpoint = fill(complex(NaN, NaN), length(p))

    for i in 1:N
        q = all(isfinite, path[i]) ? trackStep(f, g, γ, path[i], λs[i], λs[i + 1]) : nothing
        if q === nothing || !all(isfinite, q) || norm(q) > escape
            path[i + 1] = nanpoint
        else
            path[i + 1] = q
        end
    end

    λs, path
end

"""
    followHomotopy(initial_saddles, Fs, N = 10000, γ = DEFAULT_γ; ρ0 = nothing)

Track all saddle points of γ g into the saddle points of f (Fs = [f, g]). Returns the paths, their
Borel singularities t = -i h(x, λ) and, along every path, a continuous branch ρ of √det Hess h
(starting from ρ0 when given), which fixes the normalisation 1/ρ of the asymptotic series.
"""
function followHomotopy(initial_saddles, Fs, N = 10000, γ = DEFAULT_γ; ρ0 = nothing, escape = 1e6)
    f, g = Fs[1], Fs[2]
    n = length(initial_saddles)
    paths = Vector{Vector{Vector{ComplexF64}}}(undef, n)
    pathsBorel = Vector{Vector{ComplexF64}}(undef, n)
    ρs = Vector{Vector{ComplexF64}}(undef, n)

    h(p, λ) = λ * f(p) + γ * (1 - λ) * g(p)

    Threads.@threads for i in 1:n
        λs, path = homotopy(initial_saddles[i], Fs, N, γ; escape)
        paths[i] = path
        pathsBorel[i] = [all(isfinite, p) ? -im * h(p, λ) : complex(NaN, NaN) for (p, λ) in zip(path, λs)]

        # continuous branch of √det Hess h along the path
        ρ = fill(complex(NaN, NaN), length(path))
        r = sqrt(det(cHessian(p -> h(p, 0.0), path[1])))
        ρ0 === nothing || (r = abs(r - ρ0[i]) < abs(r + ρ0[i]) ? r : -r)
        ρ[1] = r
        for k in 2:length(path)
            all(isfinite, path[k]) || break
            s = sqrt(det(cHessian(p -> h(p, λs[k]), path[k])))
            ρ[k] = abs(s - ρ[k - 1]) ≤ abs(s + ρ[k - 1]) ? s : -s
        end
        ρs[i] = ρ
    end

    return paths, pathsBorel, ρs
end
