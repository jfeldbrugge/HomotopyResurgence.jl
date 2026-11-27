```@meta
CurrentModule = HomotopyResurgence
```

# Theory

HomotopyResurgence.jl evaluates the Stokes constants of oscillatory integrals

```math
\Psi(k) = \left(\frac{k}{2\pi i}\right)^{d/2} \int_{\mathbb{R}^d} e^{i k f(\boldsymbol{x})}\, \mathrm{d}\boldsymbol{x}
```

with a polynomial exponent $f$. The Stokes constants decide how the asymptotic series of the saddle points see each other in the Borel plane, and with them which saddle points contribute to the integral. Computing them from the divergent tails of the series, or from the intersections of steepest descent manifolds in $\mathbb{C}^d$, is costly. The package takes a different route: it starts from an integral whose Stokes constants are known, deforms it into the integral of interest and keeps track of the Stokes constants along the way. This page summarises the method of *Stokes Constants through Homotopy* (Aniceto, Feldbrugge and Howls, 2026).

## Stokes constants and the Borel plane

Let $\boldsymbol{x}_j$ be the non-degenerate saddle points of $f$, $\nabla f(\boldsymbol{x}_j) = 0$, with the Borel singularities $t_j = -i f(\boldsymbol{x}_j)$. Along the steepest descent manifold of $\boldsymbol{x}_j$, the integral has the asymptotic expansion

```math
\Psi^{(j)}(k) \sim e^{i k f(\boldsymbol{x}_j)} \sum_{m=0}^\infty \frac{T_m^{(j)}}{k^m}, \qquad T_0^{(j)} = \frac{1}{\sqrt{\det \mathcal{H}f(\boldsymbol{x}_j)}}.
```

The late terms of the series of $\boldsymbol{x}_j$ contain the early terms of the series of the other saddle points (Berry & Howls 1991),

```math
T_N^{(j)} \sim \frac{1}{2\pi i} \sum_l K_{jl} \sum_{m=0}^\infty \frac{(N - m - 1)!}{F_{jl}^{\,N-m}}\, T_m^{(l)}, \qquad F_{jl} = t_l - t_j,
```

with the integer Stokes constants $K_{jl}$. Two saddle points are adjacent when $K_{jl} \neq 0$. The Stokes constants depend on the branch of the square root in $T_0^{(j)}$: flipping the series of $\boldsymbol{x}_j$ replaces $K_{jl}$ by $-K_{jl}$. The package fixes the branch with [`seriesNormalisation`](@ref): the principal branch for a complex saddle point, and stationary phase, $\sqrt{\det \mathcal{H}f} = i^{\nu}\sqrt{|\det \mathcal{H}f|}$ with the Morse index $\nu$, for a real one.

## The reference integral

For the exponent $g(x) = x^m/m - x$ in one variable, the $m - 1$ saddle points $e_n = e^{2\pi i n/(m-1)}$ lie on the unit circle, and every pair is adjacent,

```math
K_{jl} = \mathrm{sgn}(l - j),
```

when the series are normalised continuously along the circle, $\rho_n = \sqrt{m-1}\, e^{-i\pi n/(m-1)}$ ([`referenceStokes1D`](@ref)). The reference integral of the package is a sum of such exponents,

```math
R(\boldsymbol{x}) = \sum_{n=1}^d c_n \left(\frac{x_n^{m_n}}{m_n} - x_n\right),
```

with generic phases $c_n$ ([`referenceSystem`](@ref), [`goldenPhases`](@ref)). Its integral factorises into one-dimensional integrals. Two of its saddle points see each other only when they differ in a single coordinate, with the Stokes constant of that coordinate ([`initialSaddles`](@ref)). The phases $c_n$ remove the degenerate critical values that would arise when two indices $m_n$ agree. The indices follow from the degrees of the partial derivatives of $f$, $m_n - 1 \geq \deg \partial f/\partial x_n$, so that the reference integral has at least as many saddle points as $f$ (Bézout's theorem).

## Homotopy continuation

The homotopy

```math
H(\boldsymbol{x}, \lambda) = \lambda f(\boldsymbol{x}) + \gamma (1 - \lambda) R(\boldsymbol{x}), \qquad \lambda \in [0, 1],
```

with a generic phase $\gamma$ (the gamma trick of total-degree homotopies), deforms the reference integral into the integral of interest. The saddle points move along tracks $\boldsymbol{x}_j(\lambda)$, which solve

```math
\frac{\mathrm{d}\boldsymbol{x}_j}{\mathrm{d}\lambda} = -\left[\mathcal{H}_{\boldsymbol{x}} H\right]^{-1} \partial_\lambda \nabla_{\boldsymbol{x}} H,
```

and end at the saddle points of $f$ or at infinity. [`homotopy`](@ref) integrates one track with a predictor–corrector scheme, an Euler step followed by Newton's method, halving the step until Newton's method converges quickly. [`followHomotopy`](@ref) tracks all saddle points and continues the branch of $\sqrt{\det \mathcal{H}_{\boldsymbol{x}} H}$ along every track, which fixes the normalisation of the series continuously in $\lambda$.

## Higher-order Stokes phenomena

The Stokes constants change only when three Borel singularities line up, with $t_p$ between $t_j$ and $t_l$: a higher-order Stokes phenomenon (Howls, Langman & Olde Daalhuis 2004). When $t_p$ crosses the directed segment from $t_j$ to $t_l$, the Stokes constant of the outer pair changes by

```math
K_{jl} \mapsto K_{jl} + \sigma K_{jp} K_{pl},
```

with $\sigma = +1$ for a crossing from left to right and $\sigma = -1$ for a crossing from right to left. The paper derives the rule from the exact remainder of the asymptotic series. [`findColinear`](@ref) finds the collinearities of three tracks, and [`higherOrderStokes`](@ref) applies the rule to all events in the order of $\lambda$. With the keyword `upto`, it returns the Stokes constants at any point of the homotopy. The signed rule carries Stokes constants of any magnitude. A bookkeeping of the adjacency alone, which toggles an edge at every event, fails as soon as a Stokes constant of magnitude two enters.

At $\lambda = 1$, [`stokesHomotopy`](@ref) converts the continuous normalisations of the tracks to the normalisation of the series of $f$, discards the saddle points that flew off to infinity and returns the Stokes constants of the saddle points of $f$.

## Rational exponents

For a rational exponent $f = P/Q$, introduce two auxiliary variables: $x_{d+1}$ replaces $1/Q$, and $x_{d+2}$ enforces $x_{d+1} Q = 1$,

```math
P_{\mathrm{ext}}(\boldsymbol{x}, x_{d+1}, x_{d+2}) = P(\boldsymbol{x})\, x_{d+1} + x_{d+2}\left(x_{d+1} Q(\boldsymbol{x}) - 1\right).
```

The polynomial $P_{\mathrm{ext}}$ has the saddle points and critical values of $f$, and the homotopy applies to it unchanged. The lift raises the degrees of the partial derivatives, so many saddle points of the reference integral fly off to infinity. On their way out, they can still cross the segments between the surviving saddle points and change their Stokes constants. The Hankel function, with $f = x + 1/x$, acquires its Stokes constant $2$ in this way, and the mirror of $\mathbb{CP}^2$, with $f = x + y + 1/(xy)$, its Stokes constants $3$.

## Bounds

For a polynomial exponent in one variable, at most one steepest descent path joins two saddle points when no third critical value lies on the segment between them, and $|K_{jl}| \leq 1$. A pole lifts this bound: the two paths of the Hankel function enclose the pole at $x = 0$. In more variables, the intersection form of the vanishing cycles bounds the Stokes constants. The constants of the mirror of $\mathbb{CP}^n$ grow as the binomial coefficients, $|K_{jl}| = \binom{n+1}{|j - l|}$.
