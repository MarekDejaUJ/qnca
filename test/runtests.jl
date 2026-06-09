using Test
using Random
include(joinpath(@__DIR__, "..", "src", "QNCA.jl"))
using .QNCA

const SCOPE = (0.0, 100.0, 0.0, 100.0)

@testset "validation invariant: d_pi(pi=1) == native CE-FDH" begin
    for n in (100, 200, 500, 1000)
        rng = MersenneTwister(42)
        X, Y = generate_reverse_L(n; rng=rng)
        for scope in (SCOPE, nothing)
            q1 = qnca(X, Y; pi=1.0, n_grid=400, B=0, scope=scope)
            ref = nca_ce_fdh_d(X, Y; scope=scope)
            @test abs(q1.d_pi - ref) < 2e-3
        end
    end
end

@testset "frontier non-decreasing" begin
    rng = MersenneTwister(1)
    X, Y = generate_reverse_L(300; rng=rng)
    y_grid = collect(range(0.0, 100.0, length=50))
    phi = qnca_frontier(X, Y, 0.95, y_grid)
    fin = filter(!isnan, phi)
    @test all(diff(fin) .>= -1e-9)
end

@testset "effect size in [0,1] and grows as pi falls" begin
    rng = MersenneTwister(2)
    X, Y = generate_reverse_L(500; rng=rng)
    d1  = qnca(X, Y; pi=1.00, n_grid=200, B=0, scope=SCOPE).d_pi
    d95 = qnca(X, Y; pi=0.95, n_grid=200, B=0, scope=SCOPE).d_pi
    d90 = qnca(X, Y; pi=0.90, n_grid=200, B=0, scope=SCOPE).d_pi
    @test 0.0 <= d1 && d90 <= 1.0
    @test d1 <= d95 + 1e-9
    @test d95 <= d90 + 1e-9
end

@testset "permutation p floor 1/(B+1)" begin
    rng = MersenneTwister(7)
    X, Y = generate_reverse_L(200; rng=rng)
    B = 199
    res = qnca(X, Y; pi=0.95, B=B, scope=SCOPE, seed=1)
    @test res.p_pi >= 1/(B+1) - 1e-12
    @test res.p_pi <= 1.0
end

@testset "supplied permutations control p-value" begin
    rng = MersenneTwister(11)
    X, Y = generate_reverse_L(50; rng=rng)
    permutations = repeat(reshape(collect(1:length(Y)), 1, :), 5, 1)
    res = qnca(X, Y; pi=0.95, B=0, scope=SCOPE, permutations=permutations)
    @test res.p_pi == 1.0
end

@testset "no signal -> small effect" begin
    rng = MersenneTwister(0)
    X = rand(rng, 500) .* 100.0
    Y = rand(rng, 500) .* 100.0
    res = qnca(X, Y; pi=1.0, n_grid=200, B=0, scope=SCOPE)
    @test res.d_pi < 0.10
end

@testset "spuriousness band shape and bounds" begin
    rng = MersenneTwister(21)
    X, Y = generate_reverse_L(80; rng=rng)
    draws = randn(MersenneTwister(3), 12, length(X), 2)
    res = spuriousness_band(X, Y; pi_grid=[1.0, 0.95], n_grid=40,
                            M=12, scope=SCOPE, normal_draws=draws)
    @test length(res.band.pi) == 2
    @test size(res.d_null) == (12, 2)
    @test all(isfinite, res.band.d_obs)
    @test all(res.band.lower .<= res.band.median .+ 1e-12)
    @test all(res.band.median .<= res.band.upper .+ 1e-12)
    @test all(res.band.p_null .>= 1 / 13 - 1e-12)
    @test all(res.band.p_null .<= 1.0)
end

@testset "spuriousness band accepts uniform draws" begin
    rng = MersenneTwister(31)
    X, Y = generate_reverse_L(60; rng=rng)
    uniforms = rand(MersenneTwister(5), 9, length(X), 2)
    res = spuriousness_band(X, Y; pi_grid=[1.0, 0.95], n_grid=35,
                            M=9, scope=SCOPE, uniform_draws=uniforms)
    @test size(res.d_null) == (9, 2)
    @test all(isfinite, res.band.d_obs)
    @test all(isfinite, res.band.lower)
    @test all(isfinite, res.band.upper)
end

@testset "spuriousness band rejects two supplied streams" begin
    rng = MersenneTwister(32)
    X, Y = generate_reverse_L(30; rng=rng)
    uniforms = rand(MersenneTwister(6), 4, length(X), 2)
    normals = randn(MersenneTwister(7), 4, length(X), 2)
    @test_throws ErrorException spuriousness_band(
        X, Y; M=4, scope=SCOPE, uniform_draws=uniforms, normal_draws=normals
    )
end

@testset "consistency probe finite summaries" begin
    rng = MersenneTwister(22)
    X, Y = generate_reverse_L(90; rng=rng)
    res = consistency_probe(X, Y; pi_pair=(1.0, 0.90),
                            sizes=[30, 60, 90], reps=8,
                            n_grid=40, scope=SCOPE, seed=4)
    @test size(res.effects) == (3, 8, 2)
    @test length(res.summary.mean_d) == 6
    @test all(isfinite, res.summary.mean_d)
    @test all(isfinite, res.summary.sd_d)
    @test all(isfinite, res.beta_drift)
    @test isfinite(res.divergence)
end
