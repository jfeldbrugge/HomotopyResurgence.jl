```@meta
CurrentModule = HomotopyResurgence
```

# Tutorial

In this tutorial, we demonstrate the main features of HomotopyResurgence.jl with four examples: the Pearcey integral, the swallowtail catastrophe, whose Stokes constants change in higher-order Stokes phenomena along the homotopy, the parabolic umbilic in two variables, and the Hankel function, a rational exponent that we treat through a polynomial lift.

The exponent enters as a Julia function of a vector `q` of complex variables. The package differentiates it with ForwardDiff, so any polynomial written with `+`, `*` and `^` works.

## The Pearcey integral

We consider the Pearcey integral

```math
\Psi(k) = \sqrt{\frac{k}{2\pi i}} \int_{\mathbb{R}} e^{i k f(q)}\, \mathrm{d}q, \qquad f(q) = q^4 + \frac{q^2}{2} + q.
```

The derivative $f'(q) = 4q^3 + q + 1$ has degree $3$, so the reference integral needs the index $m = 4$, with the $m - 1 = 3$ saddle points of $g(x) = x^4/4 - x$.

```@example pearcey
using HomotopyResurgence, CairoMakie

f(q) = q[1]^4 + q[1]^2 / 2 + q[1]
r = stokesHomotopy(f, [4]; N = 4000)
nothing # hide
```

[`stokesHomotopy`](@ref) returns a named tuple. The saddle points of $f$ that the homotopy reaches, their Borel singularities $t_j = -i f(q_j)$ and the Stokes constants are

```@example pearcey
first.(r.x)
```

```@example pearcey
r.Σ
```

```@example pearcey
r.K
```

`K[j, l]` is the Stokes constant $K_{jl}$ with which the series of saddle point $j$ sees the Borel singularity of saddle point $l$, in the convention of Berry & Howls (1991). Every pair of saddle points of the Pearcey integral is adjacent. The normalisations $T_0^{(j)} = 1/\sqrt{f''(q_j)}$ of the series, for which these constants hold, are

```@example pearcey
r.T0
```

The homotopy started from the reference integral, with the Stokes constants of [`referenceStokes1D`](@ref):

```@example pearcey
r.K0
```

The tracks of the Borel singularities, with the adjacency graph at the end of the homotopy:

```@example pearcey
fig, ax = plotHomotopy(r.pathsBorel, r.K, 1.6)
fig
```

## The swallowtail catastrophe: higher-order Stokes phenomena

For the swallowtail catastrophe,

```math
f(t) = t^5 + z t^3 + y t^2 + x t, \qquad x = y = 0.7, \quad z = -1,
```

the adjacency graph changes along the homotopy. We use the phase $\gamma = e^{i\pi/5}$ of the paper:

```@example swallowtail
using HomotopyResurgence, CairoMakie

f(q) = q[1]^5 - q[1]^3 + 0.7q[1]^2 + 0.7q[1]
r = stokesHomotopy(f, [5]; N = 4000, γ = cis(π / 5))
r.K
```

The higher-order Stokes phenomena of the homotopy are stored in `r.events`. Every event records the grid parameter `τ` (with $\lambda = (\tau - 1)/N$), the grid index `i` and the triple: the Borel singularity of `k` crosses the segment between those of `j` and `l`.

```@example swallowtail
r.events
```

[`higherOrderStokes`](@ref) replays these events. With the keyword `upto`, it returns the Stokes constants at any point of the homotopy, here before and after the first event (for the tracks of the reference integral, before the conversion to the normalisation of $f$):

```@example swallowtail
τ = r.events[1].τ
higherOrderStokes(r.K0, r.pathsBorel; upto = floor(Int, τ), events = r.events)
```

```@example swallowtail
higherOrderStokes(r.K0, r.pathsBorel; upto = ceil(Int, τ) + 1, events = r.events)
```

The adjacency graph at the start, halfway and at the end of the homotopy:

```@example swallowtail
N = length(r.pathsBorel[1]) - 1
fig = Figure(size = (1200, 400))
for (n, i) in enumerate((1, N ÷ 2 + 1, N + 1))
    K = higherOrderStokes(r.K0, r.pathsBorel; upto = i, events = r.events)
    ax = Axis(fig[1, n], title = "λ = $((i - 1) / N)", aspect = 1, limits = (-0.8, 0.8, -0.8, 0.8))
    pts = [p[i] for p in r.pathsBorel]
    for p in r.pathsBorel
        lines!(ax, real.(p[1:i]), imag.(p[1:i]), color = (:steelblue, 0.6))
    end
    for a in eachindex(pts), b in a+1:length(pts)
        K[a, b] == 0 && continue
        lines!(ax, real.(pts[[a, b]]), imag.(pts[[a, b]]), color = :red, linewidth = 2abs(K[a, b]))
    end
    scatter!(ax, real.(pts), imag.(pts), color = :black)
end
fig
```

## The parabolic umbilic: two variables

The parabolic umbilic has the exponent

```math
f(u, v) = u^2 v + v^4 + x u^2 + y v^2 + z u + w v, \qquad x = 0.1, \quad y = -0.1, \quad z = 1, \quad w = 0.
```

The partial derivatives have the degrees $2$ and $3$, so the reference integral has the indices $(3, 4)$ and $2 \cdot 3 = 6$ saddle points:

```@example parabolic
using HomotopyResurgence

f(q) = q[1]^2 * q[2] + q[2]^4 + 0.1q[1]^2 - 0.1q[2]^2 + q[1]
r = stokesHomotopy(f, [3, 4]; N = 4000)
length(r.pathsBorel), r.finite
```

Five tracks reach the saddle points of $f$, and the sixth flies off to infinity. The Stokes constants of the five saddle points:

```@example parabolic
r.K
```

## The Hankel function: a rational exponent

The Hankel function of order one,

```math
\Psi(k) = \int_0^\infty e^{i k (x + 1/x)}\, \mathrm{d}x = -\pi H_1^{(1)}(2k),
```

has a rational exponent with the saddle points $x = \pm 1$. To apply the homotopy, lift it to a polynomial exponent in three variables: $x_2$ replaces $1/x_1$, and $x_3$ enforces $x_2 x_1 = 1$,

```math
P_{\mathrm{ext}}(x_1, x_2, x_3) = (1 + x_1^2)\, x_2 + x_3\,(x_2 x_1 - 1).
```

Every partial derivative has degree $2$, so the reference integral has the indices $(3, 3, 3)$ and $8$ saddle points:

```@example hankel
using HomotopyResurgence

P(q) = (1 + q[1]^2) * q[2] + q[3] * (q[2] * q[1] - 1)
r = stokesHomotopy(P, [3, 3, 3]; N = 4000)
length(r.pathsBorel), r.finite
```

Six saddle points fly off to infinity. The two that remain reach $x_1 = \pm 1$, with the Borel singularities $t = \pm 2i$:

```@example hankel
r.Σ
```

```@example hankel
r.K
```

The Stokes constant has magnitude $2$, beyond the bound $|K_{jl}| \leq 1$ of polynomial exponents in one variable. The escaping saddle points build it: on their way to infinity, they cross the segment between the two survivors, and every crossing adds $\pm 1$.

## The normalisation of the series

The Stokes constants depend on the branch of $\sqrt{\det \mathcal{H}f}$ in the normalisation of every series. [`seriesNormalisation`](@ref) fixes it: by default (`:morse`), the principal branch for a complex saddle point and stationary phase, $\rho = i^{\nu}\sqrt{|\det \mathcal{H}f|}$ with the Morse index $\nu$, for a real one, as in `transseries` of IntegralResurgence.jl. The two conventions differ at a real maximum in two variables, with $\mathcal{H}f = -I$ and $\nu = 2$: stationary phase gives $\rho = i^2 = -1$, the principal branch $\rho = +1$.

```@example normalisation
using HomotopyResurgence
H = ComplexF64[-1 0; 0 -1]
x = ComplexF64[0, 0]
seriesNormalisation(H, x), seriesNormalisation(H, x; normalisation = :principal)
```

Pass `normalisation = :principal` to [`stokesHomotopy`](@ref) to use the principal branch for every saddle point. Flipping the branch of the series of saddle point $j$ flips the signs of row and column $j$ of `K`.
