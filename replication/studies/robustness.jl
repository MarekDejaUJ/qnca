# Robustness evaluation, studies S1 to S4 (design in README.md).
#
# Generates every sample from fixed seeds, writes the base samples, planted
# positions and noise draws for the NCA computations in R, and computes the
# CE-FDH, QNCA, order-m and screen-and-delete results.
#
#   julia -t auto replication/studies/robustness.jl [reps]

include(joinpath(@__DIR__, "..", "..", "julia", "QuantileNecessity.jl", "src", "QuantileNecessity.jl"))
using .QuantileNecessity
using Random, Printf, Base.Threads

const ROOT = normpath(joinpath(@__DIR__, ".."))
const OUT = joinpath(ROOT, "results")
const DATA = joinpath(OUT, "data")
const BASE_SEED = 20260930
const SCOPE = (0.0, 100.0, 0.0, 100.0)
const N_GRID = 200
const GRID = collect(range(0.0, 100.0, length=N_GRID))
const PI_MAIN = [0.95, 0.90]
const PI_PATH = [1.0, 0.99, 0.975, 0.95, 0.925, 0.90, 0.85, 0.80]
const ORDER_M = [19, 9]
const M_LEVELS = [1, 2, 3, 5]
const SIGMAS = [0.0, 2.0, 5.0, 10.0]
const GENERATORS = ["G1", "G2"]
const N_S1 = [100, 300, 1000]
const N_S2 = [100, 300]
const Y_BOTTLENECK = 50.0

reps = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 500

# Seed for one replicate of one cell: an explicit formula, stable across
# Julia versions.
cell_seed(study, gen, n, rep) =
    BASE_SEED + 100_000_000 * study + 10_000_000 * findfirst(==(gen), GENERATORS) +
    10_000 * n + rep

ceiling_value(gen, x) = gen == "G1" ? 15.0 + 0.85x : 10.0 * sqrt(x)

function draw_sample(rng, gen, n, noise_sd)
    X = rand(rng, n) .* 100.0
    U = rand(rng, n)
    e = randn(rng, n) .* noise_sd
    Y = (U .^ 2) .* ceiling_value.(gen, X) .+ e
    return clamp.(X, 0.0, 100.0), clamp.(Y, 0.0, 100.0)
end

function draw_planted(rng)
    mid = [(10.0 * rand(rng), 40.0 + 30.0 * rand(rng)) for _ in 1:5]
    top = [(10.0 * rand(rng), 80.0 + 20.0 * rand(rng)) for _ in 1:5]
    return mid, top
end

# ---------------------------------------------------------------- estimators

"Linear interpolation of a grid frontier at outcome level y."
function at_level(phi, y)
    j = searchsortedlast(GRID, y)
    j >= length(GRID) && return phi[end]
    t = (y - GRID[j]) / (GRID[j+1] - GRID[j])
    return (1 - t) * phi[j] + t * phi[j+1]
end

"CE-FDH exact effect and exact bottleneck r(y) = min{X : Y >= y}."
function cefdh(X, Y)
    d = nca_ce_fdh_d(X, Y; scope=SCOPE)
    hit = Y .>= Y_BOTTLENECK
    b = any(hit) ? minimum(X[hit]) : SCOPE[2]
    return d, b
end

"Raw QNCA frontiers for several tolerances from one pass over the grid."
function raw_frontiers(X, Y, pis)
    qs = [QuantileNecessity.decimal_tail_probability(p) for p in pis]
    parts = [QuantileNecessity.tail_parts(q) for q in qs]
    phis = [fill(NaN, N_GRID) for _ in pis]
    ordm = [fill(NaN, N_GRID) for _ in ORDER_M]
    ord = sortperm(Y; rev=true)
    Ys = Y[ord]; Xs = X[ord]
    for j in 1:N_GRID
        k = searchsortedlast(Ys, GRID[j]; rev=true)
        k == 0 && continue
        x = sort(Xs[1:k])
        for p in eachindex(pis)
            phis[p][j] = x[QuantileNecessity.selected_rank(k, qs[p], parts[p])]
        end
        for (a, m) in enumerate(ORDER_M)
            s = 0.0
            for i in 1:k
                w = ((k - i + 1) / k)^m - ((k - i) / k)^m
                s += x[i] * w
            end
            ordm[a][j] = s
        end
    end
    return phis, ordm
end

function fit_d_b(raw, mode)
    phi = QuantileNecessity.fit_frontier!(copy(raw), SCOPE[2], mode)
    d = QuantileNecessity.qnca_d(phi, GRID, SCOPE[1], 1e4)
    return d, at_level(phi, Y_BOTTLENECK)
end

"CE-FDH peers: cases that are the unique running minimum at their outcome."
function cefdh_candidates(X, Y)
    cand = Int[]
    for i in eachindex(X)
        lower = true
        for l in eachindex(X)
            if l != i && Y[l] >= Y[i] && X[l] <= X[i]
                lower = false
                break
            end
        end
        lower && push!(cand, i)
    end
    return cand
end

"""
Screen and delete: repeatedly delete the case whose removal raises the CE-FDH
effect most. `rule = :oracle` deletes exactly `m` cases; `rule = :onepct`
deletes while the gain is at least 1% of the current effect, at most 10 times.
Returns the final effect, bottleneck, deletions, and planted deletions.
"""
function screen_delete(X, Y, planted_ids, rule, m)
    keep = collect(eachindex(X))
    deleted = Int[]
    max_steps = rule === :oracle ? m : 10
    for _ in 1:max_steps
        Xk = X[keep]; Yk = Y[keep]
        d0 = nca_ce_fdh_d(Xk, Yk; scope=SCOPE)
        best_gain = -Inf; best = 0
        for c in cefdh_candidates(Xk, Yk)
            sel = trues(length(keep)); sel[c] = false
            g = nca_ce_fdh_d(Xk[sel], Yk[sel]; scope=SCOPE) - d0
            if g > best_gain
                best_gain = g; best = c
            end
        end
        best == 0 && break
        if rule === :onepct && !(best_gain >= 0.01 * d0 && best_gain > 0.0)
            break
        end
        push!(deleted, keep[best])
        deleteat!(keep, best)
    end
    d, b = cefdh(X[keep], Y[keep])
    return d, b, length(deleted), count(in(planted_ids), deleted)
end

"First case the screen deletes and its effect gain (for the NCA check)."
function first_screen_case(X, Y)
    d0 = nca_ce_fdh_d(X, Y; scope=SCOPE)
    best_gain = -Inf; best = 0
    for c in cefdh_candidates(X, Y)
        sel = trues(length(X)); sel[c] = false
        g = nca_ce_fdh_d(X[sel], Y[sel]; scope=SCOPE) - d0
        if g > best_gain
            best_gain = g; best = c
        end
    end
    return best, best_gain
end

# ---------------------------------------------------------------- one dataset

struct Row
    study::String; gen::String; n::Int; rep::Int; band::String; level::Float64
    method::String; pi::Float64; d::Float64; b50::Float64
    deleted::Int; planted_deleted::Int
end

function evaluate!(rows, study, gen, n, rep, band, level, X, Y, planted_ids, m_true; path=false)
    push_row(method, pi, d, b; del=-1, pdel=-1) =
        push!(rows, Row(study, gen, n, rep, band, level, method, pi, d, b, del, pdel))
    d, b = cefdh(X, Y)
    push_row("CE-FDH", 1.0, d, b)
    pis = path ? PI_PATH : vcat(1.0, PI_MAIN)
    phis, ordm = raw_frontiers(X, Y, pis)
    for (p, pi) in enumerate(pis)
        de, be = fit_d_b(phis[p], :envelope)
        push_row("QNCA envelope", pi, de, be)
        if pi in PI_MAIN
            di, bi = fit_d_b(phis[p], :isotonic)
            push_row("QNCA least squares", pi, di, bi)
        end
    end
    for (a, m) in enumerate(ORDER_M)
        dm, bm = fit_d_b(ordm[a], :envelope)
        push_row("Order-m $(m)", NaN, dm, bm)
    end
    if study == "S1"
        ds, bs, nd, np = screen_delete(X, Y, planted_ids, :oracle, m_true)
        push_row("Screen oracle m", NaN, ds, bs; del=nd, pdel=np)
        ds, bs, nd, np = screen_delete(X, Y, planted_ids, :onepct, 0)
        push_row("Screen 1% rule", NaN, ds, bs; del=nd, pdel=np)
    end
    return rows
end

# ---------------------------------------------------------------- runners

function write_csv(path, header, lines)
    open(path, "w") do io
        println(io, header)
        for l in lines
            println(io, l)
        end
    end
end

fmt(x) = isnan(x) ? "NA" : repr(x)

const CHECK_REPS = 20

function run_s1(gen, n)
    per_rep = Vector{Vector{Row}}(undef, reps)
    check_lines = [String[] for _ in 1:reps]
    clean_lines = Vector{Vector{String}}(undef, reps)
    planted_lines = Vector{Vector{String}}(undef, reps)
    @threads for rep in 1:reps
        rng = MersenneTwister(cell_seed(1, gen, n, rep))
        X, Y = draw_sample(rng, gen, n, 4.0)
        mid, top = draw_planted(rng)
        rows = Row[]
        evaluate!(rows, "S1", gen, n, rep, "clean", 0.0, X, Y, Int[], 0; path=true)
        for (band, pts) in (("mid", mid), ("top", top)), m in M_LEVELS
            Xa = vcat(X, first.(pts[1:m])); Ya = vcat(Y, last.(pts[1:m]))
            evaluate!(rows, "S1", gen, n, rep, band, float(m), Xa, Ya,
                      collect(n+1:n+m), m; path=(m in (1, 3, 5)))
            if n <= 300 && rep <= CHECK_REPS
                c, g = first_screen_case(Xa, Ya)
                push!(check_lines[rep], "$rep,$band,$m,$c,$(repr(g))")
            end
        end
        per_rep[rep] = rows
        clean_lines[rep] = ["$rep,$i,$(repr(X[i])),$(repr(Y[i]))" for i in 1:n]
        planted_lines[rep] = vcat(["$rep,mid,$j,$(repr(mid[j][1])),$(repr(mid[j][2]))" for j in 1:5],
                                  ["$rep,top,$j,$(repr(top[j][1])),$(repr(top[j][2]))" for j in 1:5])
    end
    write_csv(joinpath(DATA, "s1_clean_$(gen)_n$(n).csv"), "rep,i,X,Y", reduce(vcat, clean_lines))
    write_csv(joinpath(DATA, "s1_planted_$(gen)_n$(n).csv"), "rep,band,j,x,y", reduce(vcat, planted_lines))
    if n <= 300
        write_csv(joinpath(DATA, "s1_screencheck_$(gen)_n$(n).csv"), "rep,band,m,case,gain",
                  reduce(vcat, check_lines))
    end
    return reduce(vcat, per_rep)
end

function run_s2(gen, n)
    per_rep = Vector{Vector{Row}}(undef, reps)
    base_lines = Vector{Vector{String}}(undef, reps)
    @threads for rep in 1:reps
        rng = MersenneTwister(cell_seed(2, gen, n, rep))
        X0, Y0 = draw_sample(rng, gen, n, 0.0)
        zx = randn(rng, n); zy = randn(rng, n)
        rows = Row[]
        for s in SIGMAS
            X = clamp.(X0 .+ s .* zx, 0.0, 100.0)
            Y = clamp.(Y0 .+ s .* zy, 0.0, 100.0)
            evaluate!(rows, "S2", gen, n, rep, "noise", s, X, Y, Int[], 0)
        end
        per_rep[rep] = rows
        base_lines[rep] = ["$rep,$i,$(repr(X0[i])),$(repr(Y0[i])),$(repr(zx[i])),$(repr(zy[i]))" for i in 1:n]
    end
    write_csv(joinpath(DATA, "s2_base_$(gen)_n$(n).csv"), "rep,i,X,Y,zx,zy", reduce(vcat, base_lines))
    return reduce(vcat, per_rep)
end

mkpath(DATA)
all_rows = Row[]
t0 = time()
for gen in GENERATORS
    for n in N_S1
        append!(all_rows, run_s1(gen, n))
        @printf("S1 %s n=%d done (%.0f s)\n", gen, n, time() - t0); flush(stdout)
    end
    for n in N_S2
        append!(all_rows, run_s2(gen, n))
        @printf("S2 %s n=%d done (%.0f s)\n", gen, n, time() - t0); flush(stdout)
    end
end
lines = ["$(r.study),$(r.gen),$(r.n),$(r.rep),$(r.band),$(r.level),$(r.method),$(fmt(r.pi)),$(fmt(r.d)),$(fmt(r.b50)),$(r.deleted),$(r.planted_deleted)"
         for r in all_rows]
write_csv(joinpath(OUT, "robustness_julia.csv"),
          "study,gen,n,rep,band,level,method,pi,d,b50,deleted,planted_deleted", lines)
open(joinpath(OUT, "robustness_julia_provenance.txt"), "w") do io
    println(io, "reps = $reps")
    println(io, "julia = $(VERSION)")
    println(io, "threads = $(nthreads())")
    ver = match(r"version\s*=\s*\"([^\"]+)\"", read(joinpath(ROOT, "..", "julia", "QuantileNecessity.jl", "Project.toml"), String))
    println(io, "packages = $(ver === nothing ? "unknown" : ver.captures[1])")
    println(io, "rows = $(length(all_rows))")
    @printf(io, "seconds = %.0f\n", time() - t0)
end
@printf("wrote %d rows in %.0f s\n", length(all_rows), time() - t0)
