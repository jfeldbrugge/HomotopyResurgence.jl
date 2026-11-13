using Homotopy
using LinearAlgebra
using Test

@testset "Homotopy.jl" begin
    @testset "reference system" begin
        # x^m/m - x: K[l, j] = sign(l - j) in the continuous normalisation (checked against the
        # Gauss-Manin Stokes constants of IntegralResurgence.jl for m = 3, ..., 7)
        @test referenceStokes1D(4) == [0 -1 -1; 1 0 -1; 1 1 0]

        # 1D: all saddles mutually adjacent, ρ² = g''
        saddles, K, ρ = initialSaddles([5])
        @test length(saddles) == 4
        @test all(s -> abs(s[1]^4 - 1) < 1e-12, saddles)
        @test K == referenceStokes1D(5)
        @test all(abs(ρ[a]^2 - 4 * saddles[a][1]^3) < 1e-12 for a in eachindex(saddles))

        # 2D (3,4): prism graph, two triangles joined by three edges -> 9 edges
        _, K, _ = initialSaddles([3, 4])
        @test count(!iszero, K) ÷ 2 == 9
        @test K == -K'

        # 3D: adjacent iff the labels differ in exactly one coordinate
        saddles, K, _ = initialSaddles([3, 3, 3])
        idx(v) = findfirst(s -> all(abs.(s .- v) .< 1e-12), saddles)
        o = idx([1, 1, 1])
        @test abs(K[idx([-1, 1, 1]), o]) == 1
        @test K[idx([-1, -1, 1]), o] == 0
        @test K[idx([-1, -1, -1]), o] == 0
        @test all(K[i, i] == 0 for i in axes(K, 1))
    end

    @testset "zero_crossings" begin
        @test Homotopy.zero_crossings([1.0, -1.0]) == [1]
        @test Homotopy.zero_crossings([1.0, 0.0, -1.0]) == [2]   # one event, not two
        @test Homotopy.zero_crossings([1.0, 0.0, 1.0]) == Int[]  # touch, no crossing
        @test Homotopy.zero_crossings([1.0, 2.0, 3.0]) == Int[]
        @test Homotopy.zero_crossings([1.0, NaN, NaN]) == Int[]  # path gone to infinity
    end

    @testset "higher-order Stokes rule" begin
        # t_2 crosses the segment t_1 → t_3 from its left (Im > 0) to its right
        N = 10
        p1 = fill(0.0 + 0.0im, N + 1)
        p3 = fill(2.0 + 0.0im, N + 1)
        p2 = [1.0 + (0.5 - i / N) * im for i in 0:N]
        K = [0 1 0; 2 0 1; 0 3 0]                                   # K[l, j] = K_{l←j}
        K2 = higherOrderStokes(K, [p1, p2, p3])
        @test K2[3, 1] == 0 + K[2, 1] * K[3, 2]                      # left → right: +
        @test K2[1, 3] == 0 - K[2, 3] * K[1, 2]                      # seen from t_3: right → left
        @test all(K2[l, j] == K[l, j] for (l, j) in ((2, 1), (1, 2), (2, 3), (3, 2)))
    end

    @testset "Pearcey: Stokes constants from the homotopy" begin
        # f = q^4 + q^2/2 + q; ground truth from IntegralResurgence.stokes_constants
        # (T0 = principal 1/√f'', K[l, j] = K_{l←j}), saddles ordered by their Borel singularity
        f(p) = p[1]^4 + p[1]^2 / 2 + p[1]
        r = stokesHomotopy(f, [4]; N = 2000)
        @test length(r.finite) == 3
        @test all(x -> abs(4x[1]^3 + x[1] + 1) < 1e-8, r.x)
        Σref = [0.3125im, 0.5787580992953792 - 0.09375im, -0.5787580992953792 - 0.09375im]
        perm = [argmin(abs.(r.Σ .- σ)) for σ in Σref]
        @test allunique(perm)
        @test maximum(abs, r.Σ[perm] .- Σref) < 1e-8
        @test r.K[perm, perm] == [0 1 -1; -1 0 -1; 1 1 0]
    end
end
