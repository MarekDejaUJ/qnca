using Test
using Random
include(joinpath(@__DIR__, "..", "src", "QNCA.jl"))
using .QNCA

const SCOPE = (0.0, 100.0, 0.0, 100.0)

@testset "decimal tolerance selects the stated type-1 rank" begin
    @test qnca_rank(100, 0.95) == 5
    @test qnca_rank(20, 0.95) == 1
    @test qnca_rank(21, 0.95) == 2
    @test qnca_rank(10, 0.90) == 1
    @test qnca_rank(11, 0.90) == 2
    @test qnca_rank(100, 1.0) == 1
    @test qnca_rank(100, 19 // 20) == 5
    @test quantile_type1_pi(collect(1.0:100.0), 0.95) == 5.0
    for k in 1:500
        @test qnca_rank(k, 0.95) == max(1, ceil(Int, k // 20))
        @test qnca_rank(k, 0.90) == max(1, ceil(Int, k // 10))
    end
    @test_throws ArgumentError qnca_rank(100, 0.0)
    @test_throws ArgumentError qnca_rank(100, 1.01)
end

@testset "orientation sentinel distinguishes X from Y" begin
    X = [0.10, 0.20, 0.40, 0.80, 0.90]
    Y = [0.10, 0.20, 0.60, 0.70, 0.95]
    condition_required = only(qnca_frontier(X, Y, 1.0, [0.60]))
    swapped_required = only(qnca_frontier(Y, X, 1.0, [0.60]))
    @test condition_required == 0.40
    @test swapped_required == 0.70
    @test condition_required != swapped_required
end

@testset "frontier uses the corrected rank at a decimal boundary" begin
    X = collect(1.0:100.0)
    Y = collect(1.0:100.0)
    phi = qnca_frontier(X, Y, 0.95, [1.0])
    @test phi == [5.0]
end

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

@testset "generated permutation stream is reproducible" begin
    rng = MersenneTwister(8)
    X, Y = generate_reverse_L(120; rng=rng)
    a = qnca(X, Y; pi=0.95, B=79, scope=SCOPE, seed=20260830)
    b = qnca(X, Y; pi=0.95, B=79, scope=SCOPE, seed=20260830)
    @test a.p_pi == b.p_pi
    @test a.d_pi == b.d_pi
end

@testset "invalid inputs and scopes fail explicitly" begin
    @test_throws ArgumentError qnca([1.0, 2.0], [1.0]; B=0)
    @test_throws ArgumentError qnca([1.0, NaN], [1.0, 2.0]; B=0)
    @test_throws ArgumentError qnca([1.0, 1.0], [1.0, 2.0]; B=0)
    @test_throws ArgumentError qnca([1.0, 2.0], [1.0, 2.0]; B=0,
                                    scope=(0.0, 1.5, 0.0, 2.0))
    @test_throws ArgumentError qnca([1.0, 2.0], [1.0, 2.0]; B=0, n_grid=1)
end

@testset "zero-valued cohorts remain valid observations" begin
    X = [0.0, 0.20, 0.40, 0.80]
    Y = [0.0, 0.25, 0.50, 0.75]
    result = qnca(X, Y; pi=0.95, n_grid=101, B=0,
                  scope=(0.0, 1.0, 0.0, 1.0))
    @test isfinite(result.d_pi)
    @test 0.0 <= result.d_pi <= 1.0
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

@testset "envelope leaves the strict member unchanged" begin
    rng = MersenneTwister(401)
    for n in (60, 300)
        X, Y = generate_reverse_L(n; rng=rng)
        a = qnca(X, Y; pi=1.0, n_grid=120, B=0, scope=SCOPE, monotone=:envelope)
        b = qnca(X, Y; pi=1.0, n_grid=120, B=0, scope=SCOPE, monotone=:isotonic)
        @test a.x_required == b.x_required
        @test a.d_pi == b.d_pi
    end
    @test monotone_envelope([3.0, 1.0, 4.0, 2.0, 5.0]) == [3.0, 3.0, 4.0, 4.0, 5.0]
    @test_throws ArgumentError qnca([1.0, 2.0], [1.0, 2.0]; B=0, monotone=:median)
end

@testset "tolerance path is monotone in pi" begin
    rng = MersenneTwister(402)
    grid = [1.0, 0.99, 0.975, 0.95, 0.925, 0.90, 0.85, 0.80, 0.50]
    for rep in 1:20, mode in (:envelope, :isotonic)
        X, Y = generate_reverse_L(150; rng=rng)
        d = [qnca(X, Y; pi=p, n_grid=60, B=0, scope=SCOPE, monotone=mode).d_pi for p in grid]
        @test all(diff(d) .>= -1e-12)
    end
end

@testset "added cases cannot push the raw ordinate below the stated order statistic" begin
    rng = MersenneTwister(403)
    grid = collect(range(0.0, 100.0, length=40))
    for rep in 1:40
        X, Y = generate_reverse_L(200; rng=rng)
        c = rand(rng, 1:4)
        Xa = vcat(X, 20 .* rand(rng, c)); Ya = vcat(Y, 40 .+ 60 .* rand(rng, c))
        for p in (0.95, 0.90)
            q = QNCA.decimal_tail_probability(p)
            _, k, h = QNCA.raw_frontier(X, Y, q, grid)
            raw_a, _, _ = QNCA.raw_frontier(Xa, Ya, q, grid)
            env_a = monotone_envelope(raw_a[.!isnan.(raw_a)])
            bound = fill(-Inf, length(grid))
            for j in eachindex(grid)
                k[j] == 0 && continue
                xs = sort(X[Y .>= grid[j]])
                if h[j] > c
                    @test raw_a[j] >= xs[h[j] - c]
                    bound[j] = xs[h[j] - c]
                end
            end
            running = accumulate(max, bound)
            for j in eachindex(env_a)
                @test env_a[j] >= running[j]
            end
        end
    end
end

@testset "resolution table agrees with the fitted frontier" begin
    rng = MersenneTwister(404)
    X, Y = generate_reverse_L(250; rng=rng)
    Xa = vcat(X, 2.0); Ya = vcat(Y, 95.0)
    for mode in (:envelope, :isotonic)
        res = qnca_resolution(Xa, Ya; pi=0.95, n_grid=80, scope=SCOPE, monotone=mode)
        fit = qnca(Xa, Ya; pi=0.95, n_grid=80, B=0, scope=SCOPE, monotone=mode)
        @test res.fitted == fit.x_required
        ok = res.k .> 0
        @test all(res.rank[ok] .== [qnca_rank(kk, 0.95) for kk in res.k[ok]])
        @test all(res.resistance[ok] .== res.rank[ok] .- 1)
        if mode === :envelope
            @test all(res.below[ok .& .!res.inherited] .<= res.rank[ok .& .!res.inherited] .- 1)
            @test any(res.inherited)
        end
    end
end

@testset "deletion screen across the tolerance grid" begin
    rng = MersenneTwister(405)
    X, Y = generate_reverse_L(600; rng=rng)
    Xa = vcat(X, 1.0); Ya = vcat(Y, 50.0)
    out = qnca_outliers(Xa, Ya; pi_grid=[1.0, 0.95], scope=SCOPE, n_grid=60)
    @test length(out.cases) == 601
    @test out.cases[1] == [601]
    @test out.flagged[1]
    @test out.dif_abs[1, 1] > 10 * abs(out.dif_abs[1, 2])
    top = qnca_outliers(vcat(X[1:150], 1.0), vcat(Y[1:150], 95.0);
                        pi_grid=[1.0, 0.95], scope=SCOPE, n_grid=60)
    @test top.cases[1] == [151]
    @test top.dif_abs[1, 1] > top.dif_abs[1, 2] > 0.0
    joint = qnca_outliers(Xa[1:200], Ya[1:200]; pi_grid=[1.0, 0.95], k=2, scope=SCOPE,
                          n_grid=60, max_candidates=6)
    @test length(joint.cases) == 15
    @test all(length.(joint.cases) .== 2)
end

@testset "breakdown curve under added efficient cases" begin
    rng = MersenneTwister(406)
    X, Y = generate_reverse_L(300; rng=rng)
    res = qnca_sensitivity(X, Y; points=(3.0, 60.0), counts=0:6, scope=SCOPE, n_grid=80)
    @test all(res.retention[1, :] .== 1.0)
    @test all(diff(res.d; dims=1) .<= 1e-12)
    @test res.retention[2, 1] < res.retention[2, 2]
    @test_throws ArgumentError qnca_sensitivity(X, Y; points=(150.0, 60.0), scope=SCOPE)
end
