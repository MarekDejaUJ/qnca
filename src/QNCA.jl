"""
    QNCA

Quantile Necessary Condition Analysis: a tolerance-parameterised generalisation
of Dul's Necessary Condition Analysis. At each outcome target the frontier is a
low type-1 quantile of the condition among the cases attaining the target; the
deterministic NCA CE-FDH ceiling is the `pi == 1` member of these
tolerance-indexed frontiers.

The raw quantile frontier is made non-decreasing by its monotone envelope, the
running maximum over outcome targets (the default), or by least-squares isotonic
projection (`monotone = :isotonic`, the 0.3 behaviour). The implementation uses
the same inverse-empirical-CDF quantile, monotone step, and carry-forward of the
frontier to the scope ceiling as the R and Python siblings. Configured decimal
tolerances are converted to exact rational ranks before frontier evaluation.
"""
module QNCA

using Random, Statistics, Printf, Base.Threads

export generate_reverse_L, isotonic_increasing, monotone_envelope,
       qnca_frontier, qnca_d, nca_ce_fdh_d, qnca, qnca_rank, quantile_type1_pi,
       qnca_resolution, qnca_outliers, qnca_sensitivity,
       spuriousness_band, consistency_probe

"Generate a reverse-L necessity scatter (condition X, outcome Y), clamped to [0,100]."
function generate_reverse_L(n; ci=15.0, cs=0.85, k=2.0, noise=4.0, rng=Random.GLOBAL_RNG)
    X = rand(rng, n) .* 100.0
    ceiling = ci .+ cs .* X
    eff = rand(rng, n) .^ k
    Y = eff .* ceiling .+ randn(rng, n) .* noise
    return clamp.(X, 0.0, 100.0), clamp.(Y, 0.0, 100.0)
end

"Inverse empirical CDF quantile (R type 1): smallest order statistic x with F(x) >= prob."
function quantile_type1(values::AbstractVector{<:Real}, prob::Float64)
    isempty(values) && return NaN
    x = sort(collect(float.(values)))
    prob <= 0.0 && return x[1]
    prob >= 1.0 && return x[end]
    idx = clamp(ceil(Int, length(x) * prob), 1, length(x))
    return x[idx]
end

function decimal_rational(value::Real)
    value isa Rational &&
        return BigInt(Base.numerator(value)) // BigInt(Base.denominator(value))
    s = string(value)
    parsed = match(r"^([+-]?)([0-9]+)(?:\.([0-9]*))?(?:[eE]([+-]?[0-9]+))?$", s)
    parsed === nothing && throw(ArgumentError("value must have a decimal representation"))
    sign_part, integer_part, fraction_part, exponent_part = parsed.captures
    fraction = fraction_part === nothing ? "" : fraction_part
    exponent = exponent_part === nothing ? 0 : parse(Int, exponent_part)
    decimal_numerator = parse(BigInt, integer_part * fraction)
    sign_part == "-" && (decimal_numerator = -decimal_numerator)
    scale = length(fraction) - exponent
    return scale >= 0 ? decimal_numerator // big(10)^scale :
                        (decimal_numerator * big(10)^(-scale)) // big(1)
end

"Return `1-pi` as an exact rational from the configured decimal representation."
function decimal_tail_probability(pi::Real)
    isfinite(pi) || throw(ArgumentError("pi must be finite"))
    p = decimal_rational(pi)
    (0 < p <= 1) ||
        throw(ArgumentError("pi must be in (0, 1]"))
    return 1 - p
end

"Selected type-1 order-statistic rank `max(1, ceil(k*(1-pi)))`."
function qnca_rank(k::Integer, pi::Real)
    k > 0 || throw(ArgumentError("conditioning-set size must be positive"))
    q = decimal_tail_probability(pi)
    return clamp(ceil(Int, k * q), 1, k)
end

"Type-1 quantile parameterised by QNCA tolerance `pi`."
function quantile_type1_pi(values::AbstractVector{<:Real}, pi::Real)
    isempty(values) && return NaN
    x = sort(collect(float.(values)))
    return x[qnca_rank(length(x), pi)]
end

function validate_xy(X, Y)
    length(X) == length(Y) || throw(ArgumentError("X and Y must have equal lengths"))
    length(X) >= 2 || throw(ArgumentError("X and Y must contain at least two observations"))
    all(isfinite, X) || throw(ArgumentError("X must contain only finite values"))
    all(isfinite, Y) || throw(ArgumentError("Y must contain only finite values"))
    minimum(X) < maximum(X) || throw(ArgumentError("X must be non-degenerate"))
    minimum(Y) < maximum(Y) || throw(ArgumentError("Y must be non-degenerate"))
    return nothing
end

function validated_scope(X, Y, scope)
    bounds = scope === nothing ?
        (minimum(X), maximum(X), minimum(Y), maximum(Y)) : Tuple(float.(scope))
    length(bounds) == 4 || throw(ArgumentError("scope must contain x_min, x_max, y_min, y_max"))
    x_min, x_max, y_min, y_max = bounds
    all(isfinite, bounds) || throw(ArgumentError("scope bounds must be finite"))
    x_min < x_max || throw(ArgumentError("scope requires x_min < x_max"))
    y_min < y_max || throw(ArgumentError("scope requires y_min < y_max"))
    minimum(X) >= x_min && maximum(X) <= x_max ||
        throw(ArgumentError("X observations must lie inside scope"))
    minimum(Y) >= y_min && maximum(Y) <= y_max ||
        throw(ArgumentError("Y observations must lie inside scope"))
    return (x_min, x_max, y_min, y_max)
end

"Pool-adjacent-violators isotonic regression (non-decreasing, unit weights)."
function isotonic_increasing(y::AbstractVector{<:Real})
    vals = Float64[]; wts = Float64[]
    for yi in y
        push!(vals, float(yi)); push!(wts, 1.0)
        while length(vals) > 1 && vals[end-1] > vals[end]
            v2 = pop!(vals); w2 = pop!(wts)
            v1 = pop!(vals); w1 = pop!(wts)
            push!(vals, (v1*w1 + v2*w2)/(w1+w2)); push!(wts, w1+w2)
        end
    end
    out = similar(float.(y))
    idx = 1
    for b in eachindex(vals)
        for _ in 1:Int(wts[b])
            out[idx] = vals[b]; idx += 1
        end
    end
    return out
end

"""
    monotone_envelope(v)

Least non-decreasing majorant of `v`: the running maximum. Applied to the raw
quantile frontier, a requirement established at a lower outcome target is
carried to every higher target, so a sparse high target cannot lower it.
"""
function monotone_envelope(v::AbstractVector{<:Real})
    out = collect(float.(v))
    @inbounds for i in 2:length(out)
        if out[i] < out[i-1]
            out[i] = out[i-1]
        end
    end
    return out
end

function monotone_mode(monotone)
    mode = Symbol(monotone)
    mode in (:envelope, :isotonic) ||
        throw(ArgumentError("monotone must be :envelope or :isotonic"))
    return mode
end

function tail_parts(q::Rational{BigInt})
    a = numerator(q); d = denominator(q)
    (a <= typemax(Int64) && d <= typemax(Int64)) || return nothing
    return (Int128(a), Int128(d))
end

function selected_rank(k::Int, q::Rational{BigInt}, parts)
    r = parts === nothing ? cld(big(k) * numerator(q), denominator(q)) :
                            cld(Int128(k) * parts[1], parts[2])
    return clamp(Int(r), 1, k)
end

"Raw type-1 quantile frontier with conditioning-set sizes and selected ranks."
function raw_frontier(X, Y, q::Rational{BigInt}, y_grid)
    m = length(y_grid)
    phi = fill(NaN, m); k = zeros(Int, m); h = zeros(Int, m)
    parts = tail_parts(q)
    @inbounds for j in 1:m
        idx = findall(>=(y_grid[j]), Y)
        if !isempty(idx)
            x = sort(X[idx])
            k[j] = length(x)
            h[j] = selected_rank(k[j], q, parts)
            phi[j] = x[h[j]]
        end
    end
    return phi, k, h
end

function fit_frontier!(phi, x_max, mode::Symbol)
    m = length(phi)
    fin = findall(!isnan, phi)
    if length(fin) >= 2
        phi[fin] = mode === :isotonic ? isotonic_increasing(phi[fin]) :
                                        monotone_envelope(phi[fin])
    end
    if !isempty(fin)
        last = fin[end]
        if last < m
            phi[last+1:m] .= x_max === nothing ? phi[last] : x_max
        end
    end
    return phi
end

"""
    qnca_frontier(X, Y, pi, y_grid; x_max=nothing, monotone=:envelope)

Quantile necessity frontier phi_pi(y), made non-decreasing and carried forward
over trailing (high-outcome) NaN levels. `monotone = :envelope` takes the running
maximum of the raw frontier; `monotone = :isotonic` applies least-squares
isotonic projection. Under a fixed scope, pass `x_max` so the whole high-outcome
band is counted as empty, matching CE-FDH.
"""
function qnca_frontier(X, Y, pi, y_grid; x_max=nothing, monotone=:envelope)
    Xf = collect(float.(X)); Yf = collect(float.(Y)); grid = collect(float.(y_grid))
    validate_xy(Xf, Yf)
    isempty(grid) && throw(ArgumentError("y_grid must not be empty"))
    all(isfinite, grid) || throw(ArgumentError("y_grid must contain only finite values"))
    issorted(grid) || throw(ArgumentError("y_grid must be non-decreasing"))
    q = decimal_tail_probability(pi)
    return qnca_frontier_validated(Xf, Yf, q, grid; x_max=x_max,
                                   monotone=monotone_mode(monotone))
end

function qnca_frontier_validated(X, Y, q::Rational{BigInt}, y_grid; x_max=nothing,
                                 monotone::Symbol=:envelope)
    phi, _, _ = raw_frontier(X, Y, q, y_grid)
    return fit_frontier!(phi, x_max, monotone)
end

"Effect size d_pi: trapezoidal empty-zone area left of the frontier / scope."
function qnca_d(phi, y_grid, x_min, scope)
    ok = findall(!isnan, phi)
    (length(ok) < 2 || scope <= 0.0) && return 0.0
    area = 0.0
    @inbounds for t in 1:length(ok)-1
        a = ok[t]; b = ok[t+1]
        wa = max(0.0, phi[a] - x_min)
        wb = max(0.0, phi[b] - x_min)
        area += (y_grid[b] - y_grid[a]) * (wa + wb) / 2
    end
    return clamp(area / scope, 0.0, 1.0)
end

function validate_permutations(permutations, n)
    perm = permutations isa AbstractVector ? reshape(Int.(permutations), 1, :) : Int.(permutations)
    ndims(perm) == 2 || error("permutations must be a B x n matrix")
    size(perm, 2) == n || error("permutations must be a B x n matrix")
    expected = collect(1:n)
    for b in 1:size(perm, 1)
        sort(vec(perm[b, :])) == expected ||
            error("each permutation row must contain 1:n exactly once")
    end
    return perm
end

function normcdf(z::Real)
    x = float(z)
    x < 0.0 && return 1.0 - normcdf(-x)
    t = 1.0 / (1.0 + 0.2316419 * x)
    poly = (((((1.330274429 * t - 1.821255978) * t) + 1.781477937) * t -
             0.356563782) * t + 0.319381530) * t
    density = exp(-0.5 * x * x) / sqrt(2.0 * pi)
    return 1.0 - density * poly
end

function norminv(p::Real)
    q = clamp(float(p), 1e-12, 1.0 - 1e-12)
    a = [-3.969683028665376e1, 2.209460984245205e2,
         -2.759285104469687e2, 1.383577518672690e2,
         -3.066479806614716e1, 2.506628277459239e0]
    b = [-5.447609879822406e1, 1.615858368580409e2,
         -1.556989798598866e2, 6.680131188771972e1,
         -1.328068155288572e1]
    c = [-7.784894002430293e-3, -3.223964580411365e-1,
         -2.400758277161838e0, -2.549732539343734e0,
          4.374664141464968e0, 2.938163982698783e0]
    d = [7.784695709041462e-3, 3.224671290700398e-1,
         2.445134137142996e0, 3.754408661907416e0]
    plow = 0.02425
    phigh = 1.0 - plow
    if q < plow
        r = sqrt(-2.0 * log(q))
        return (((((c[1] * r + c[2]) * r + c[3]) * r + c[4]) * r + c[5]) * r + c[6]) /
               ((((d[1] * r + d[2]) * r + d[3]) * r + d[4]) * r + 1.0)
    elseif q <= phigh
        r = q - 0.5
        s = r * r
        return (((((a[1] * s + a[2]) * s + a[3]) * s + a[4]) * s + a[5]) * s + a[6]) * r /
               (((((b[1] * s + b[2]) * s + b[3]) * s + b[4]) * s + b[5]) * s + 1.0)
    else
        r = sqrt(-2.0 * log(1.0 - q))
        return -(((((c[1] * r + c[2]) * r + c[3]) * r + c[4]) * r + c[5]) * r + c[6]) /
                ((((d[1] * r + d[2]) * r + d[3]) * r + d[4]) * r + 1.0)
    end
end

function average_ranks(values)
    x = collect(float.(values))
    ord = sortperm(x; alg=MergeSort)
    ranks = similar(x)
    i = 1
    while i <= length(x)
        j = i + 1
        while j <= length(x) && x[ord[j]] == x[ord[i]]
            j += 1
        end
        ranks[ord[i:j-1]] .= (i + j - 1) / 2.0
        i = j
    end
    return ranks
end

function rank_normal_correlation(X, Y)
    n = length(X)
    ax = norminv.(average_ranks(X) ./ (n + 1.0))
    ay = norminv.(average_ranks(Y) ./ (n + 1.0))
    (std(ax) == 0.0 || std(ay) == 0.0) && return 0.0
    r = cor(ax, ay)
    isfinite(r) || return 0.0
    return clamp(r, -0.999999, 0.999999)
end

function empirical_quantile_values(values, probs)
    x = sort(collect(float.(values)))
    out = Vector{Float64}(undef, length(probs))
    @inbounds for i in eachindex(probs)
        idx = clamp(ceil(Int, length(x) * clamp(float(probs[i]), 0.0, 1.0)), 1, length(x))
        out[i] = x[idx]
    end
    return out
end

function quantile_type7(values, prob)
    x = sort(collect(float.(values)))
    isempty(x) && return NaN
    (length(x) == 1 || prob <= 0.0) && return x[1]
    prob >= 1.0 && return x[end]
    h = 1.0 + (length(x) - 1.0) * prob
    j = floor(Int, h)
    gamma = h - j
    j >= length(x) && return x[end]
    return (1.0 - gamma) * x[j] + gamma * x[j + 1]
end

function scope_tuple(X, Y, scope)
    return validated_scope(X, Y, scope)
end

function d_for_pi(X, Y, pi, scope, n_grid, monotone=:envelope)
    return qnca(X, Y; pi=pi, n_grid=n_grid, B=0, scope=scope, monotone=monotone).d_pi
end

"""
    nca_ce_fdh_d(X, Y; scope=nothing)

Native CE-FDH effect size, computed the way Dul's NCA defines it: build the
upper-left free-disposal hull from the data and integrate the empty zone exactly.
Equals `qnca`'s d at pi = 1 up to grid discretization. This is the validation
invariant.
"""
function nca_ce_fdh_d(X, Y; scope=nothing)
    X = collect(float.(X)); Y = collect(float.(Y))
    validate_xy(X, Y)
    x_min, x_max, y_min, y_max = validated_scope(X, Y, scope)
    scope_area = (x_max - x_min) * (y_max - y_min)
    scope_area <= 0.0 && return 0.0

    ord = sortperm(Y; rev=true)
    Ys = Y[ord]; Xs = X[ord]
    run_min = accumulate(min, Xs)          # run_min[j] = min(Xs[1:j])
    clamp_y(v) = min(max(v, y_min), y_max)

    area = 0.0
    top_lo = clamp_y(Ys[1])
    if y_max > top_lo
        area += (y_max - top_lo) * max(0.0, x_max - x_min)
    end
    for k in 2:length(Ys)
        hi = clamp_y(Ys[k-1]); lo = clamp_y(Ys[k])
        if hi > lo
            area += (hi - lo) * max(0.0, run_min[k-1] - x_min)
        end
    end
    bot_hi = clamp_y(Ys[end])
    if bot_hi > y_min
        area += (bot_hi - y_min) * max(0.0, run_min[end] - x_min)
    end
    return clamp(area / scope_area, 0.0, 1.0)
end

"""
    qnca(X, Y; pi=1.0, n_grid=50, B=1999, scope=nothing, seed=nothing,
         permutations=nothing)

Fit QNCA at a single tolerance. The permutation test is multithreaded; set
`B=0` to skip it. Supplying a one-based B x n permutation matrix uses those rows
instead of drawing random permutations. Returns a named tuple with pi, d_pi, p_pi,
the outcome grid and the required-X frontier.
"""
function qnca(X, Y; pi=1.0, n_grid=50, B=1999, scope=nothing, seed=nothing,
              permutations=nothing, monotone=:envelope)
    X = collect(float.(X)); Y = collect(float.(Y))
    validate_xy(X, Y)
    n_grid >= 2 || throw(ArgumentError("n_grid must be at least 2"))
    B >= 0 || throw(ArgumentError("B must be non-negative"))
    q = decimal_tail_probability(pi)
    mode = monotone_mode(monotone)
    x_min, x_max, y_min, y_max = validated_scope(X, Y, scope)
    scope_area = (x_max - x_min) * (y_max - y_min)
    y_grid = collect(range(y_min, y_max, length=n_grid))
    phi_obs = qnca_frontier_validated(X, Y, q, y_grid; x_max=x_max, monotone=mode)
    d_obs = qnca_d(phi_obs, y_grid, x_min, scope_area)

    p_pi = NaN
    if permutations !== nothing
        perm = validate_permutations(permutations, length(Y))
        hits = zeros(UInt8, size(perm, 1))
        @threads for b in 1:size(perm, 1)
            Yp = Y[vec(perm[b, :])]
            d_p = qnca_d(qnca_frontier_validated(X, Yp, q, y_grid; x_max=x_max,
                                                 monotone=mode),
                          y_grid, x_min, scope_area)
            if d_p >= d_obs
                hits[b] = 1
            end
        end
        p_pi = (1 + sum(hits)) / (size(perm, 1) + 1)
    elseif B > 0
        seed0 = seed === nothing ? 42 : seed
        rng = MersenneTwister(seed0)
        perm = Matrix{Int}(undef, B, length(Y))
        for b in 1:B
            perm[b, :] = randperm(rng, length(Y))
        end
        hits = zeros(UInt8, B)
        @threads for b in 1:B
            Yp = Y[vec(perm[b, :])]
            d_p = qnca_d(qnca_frontier_validated(X, Yp, q, y_grid; x_max=x_max,
                                                 monotone=mode),
                          y_grid, x_min, scope_area)
            hits[b] = d_p >= d_obs ? 1 : 0
        end
        p_pi = (1 + sum(hits)) / (B + 1)
    end
    return (pi=pi, d_pi=d_obs, p_pi=p_pi, y_grid=y_grid, x_required=phi_obs,
            scope=scope_area, monotone=mode)
end

"""
    spuriousness_band(X, Y; pi_grid=[1.0, 0.95, 0.90], n_grid=50, M=999,
                      scope=nothing, seed=nothing, normal_draws=nothing,
                      uniform_draws=nothing)

Proposed pi-resolved spuriousness band. The diagnostic compares the observed
`d_pi` curve with a Gaussian-copula null that preserves the empirical marginals
and observed normal-score rank correlation. Supplying `uniform_draws` maps those
uniforms through the empirical marginals directly for parity runs. It does not
modify the QNCA estimator.
"""
function spuriousness_band(X, Y; pi_grid=[1.0, 0.95, 0.90], n_grid=50, M=999,
                           scope=nothing, seed=nothing, normal_draws=nothing,
                           uniform_draws=nothing, monotone=:envelope)
    X = collect(float.(X)); Y = collect(float.(Y))
    validate_xy(X, Y)
    n_grid >= 2 || throw(ArgumentError("n_grid must be at least 2"))
    pi_vals = collect(float.(pi_grid))
    isempty(pi_vals) && throw(ArgumentError("pi_grid must not be empty"))
    decimal_tail_probability.(pi_vals)
    mode = monotone_mode(monotone)
    n = length(X)
    fixed_scope = scope_tuple(X, Y, scope)
    if normal_draws !== nothing && uniform_draws !== nothing
        error("supply at most one of normal_draws and uniform_draws")
    end
    if uniform_draws !== nothing
        uniforms = float.(uniform_draws)
        ndims(uniforms) == 3 || error("uniform_draws must have dimensions M x n x 2")
        size(uniforms, 2) == n || error("uniform_draws must have dimensions M x n x 2")
        size(uniforms, 3) == 2 || error("uniform_draws must have dimensions M x n x 2")
        M = size(uniforms, 1)
        draws = nothing
    elseif normal_draws === nothing
        rng = MersenneTwister(seed === nothing ? 42 : seed)
        draws = randn(rng, M, n, 2)
        uniforms = nothing
    else
        draws = float.(normal_draws)
        ndims(draws) == 3 || error("normal_draws must have dimensions M x n x 2")
        size(draws, 2) == n || error("normal_draws must have dimensions M x n x 2")
        size(draws, 3) == 2 || error("normal_draws must have dimensions M x n x 2")
        M = size(draws, 1)
        uniforms = nothing
    end
    M > 0 || error("M must be positive")

    d_obs = [d_for_pi(X, Y, pi, fixed_scope, n_grid, mode) for pi in pi_vals]
    d_null = Matrix{Float64}(undef, M, length(pi_vals))
    r = rank_normal_correlation(X, Y)
    scale = sqrt(max(0.0, 1.0 - r * r))
    @threads for j in 1:M
        if uniforms !== nothing
            Xj = empirical_quantile_values(X, vec(uniforms[j, :, 1]))
            Yj = empirical_quantile_values(Y, vec(uniforms[j, :, 2]))
        else
            zx = vec(draws[j, :, 1])
            zy = r .* zx .+ scale .* vec(draws[j, :, 2])
            Xj = empirical_quantile_values(X, normcdf.(zx))
            Yj = empirical_quantile_values(Y, normcdf.(zy))
        end
        for k in eachindex(pi_vals)
            d_null[j, k] = d_for_pi(Xj, Yj, pi_vals[k], fixed_scope, n_grid, mode)
        end
    end

    lower = [quantile_type7(d_null[:, k], 0.025) for k in eachindex(pi_vals)]
    median = [quantile_type7(d_null[:, k], 0.500) for k in eachindex(pi_vals)]
    upper = [quantile_type7(d_null[:, k], 0.975) for k in eachindex(pi_vals)]
    excess = d_obs .- upper
    p_null = [(1.0 + count(d_null[:, k] .>= d_obs[k])) / (M + 1.0) for k in eachindex(pi_vals)]
    band = (pi=pi_vals, d_obs=d_obs, lower=lower, median=median,
            upper=upper, excess=excess, p_null=p_null)
    return (band=band, d_null=d_null, rank_correlation=r)
end

"""
    consistency_probe(X, Y; pi_pair=(1.0, 0.90), sizes=nothing, reps=100,
                      n_grid=50, scope=nothing, seed=nothing, subsamples=nothing)

Proposed extreme-vs-consistent consistency probe. The diagnostic tracks `d_pi`
across increasing subsample sizes for `pi = 1` and a lower tolerance. It is
descriptive and can be confounded by hard support bounds or marginal skewness.
"""
function consistency_probe(X, Y; pi_pair=(1.0, 0.90), sizes=nothing, reps=100,
                           n_grid=50, scope=nothing, seed=nothing, subsamples=nothing,
                           monotone=:envelope)
    X = collect(float.(X)); Y = collect(float.(Y))
    n = length(X)
    pi_vals = collect(float.(pi_pair))
    length(pi_vals) == 2 || error("pi_pair must contain exactly two tolerances")
    if sizes === nothing
        raw = [max(2, n ÷ 4), max(2, n ÷ 2), max(2, 3n ÷ 4), n]
        size_vals = unique(clamp.(raw, 2, n))
    else
        size_vals = collect(Int.(sizes))
    end
    length(size_vals) >= 2 || error("at least two subsample sizes are required")
    all((size_vals .>= 2) .& (size_vals .<= n)) ||
        error("subsample sizes must be between 2 and n")
    reps > 0 || error("reps must be positive")
    if subsamples !== nothing
        length(subsamples) == length(size_vals) ||
            error("subsamples must contain one matrix per subsample size")
    end

    fixed_scope = scope_tuple(X, Y, scope)
    effects = Array{Float64}(undef, length(size_vals), reps, length(pi_vals))
    rng = MersenneTwister(seed === nothing ? 42 : seed)
    for k in eachindex(size_vals)
        m = size_vals[k]
        if subsamples === nothing
            rows = Matrix{Int}(undef, reps, m)
            for rr in 1:reps
                rows[rr, :] = randperm(rng, n)[1:m]
            end
        else
            rows = Int.(subsamples[k])
            size(rows) == (reps, m) || error("each subsample matrix must be reps x size")
            (minimum(rows) >= 1 && maximum(rows) <= n) ||
                error("subsample indices must be one-based and within 1:n")
        end
        for rr in 1:reps
            idx = vec(rows[rr, :])
            for p in eachindex(pi_vals)
                effects[k, rr, p] = d_for_pi(X[idx], Y[idx], pi_vals[p], fixed_scope,
                                             n_grid, monotone)
            end
        end
    end

    mean_d = Matrix{Float64}(undef, length(size_vals), length(pi_vals))
    sd_d = Matrix{Float64}(undef, length(size_vals), length(pi_vals))
    for k in eachindex(size_vals), p in eachindex(pi_vals)
        vals = effects[k, :, p]
        mean_d[k, p] = mean(vals)
        sd_d[k, p] = reps > 1 ? std(vals) : 0.0
    end

    g = 1.0 ./ float.(size_vals)
    denom = sum((g .- mean(g)).^2)
    beta = fill(NaN, length(pi_vals))
    if denom > 0.0
        for p in eachindex(pi_vals)
            ybar = mean_d[:, p]
            beta[p] = sum((g .- mean(g)) .* (ybar .- mean(ybar))) / denom
        end
    end
    divergence = all(isfinite, beta) ? abs(beta[1]) - abs(beta[2]) : NaN
    summary = (size=repeat(size_vals, outer=length(pi_vals)),
               pi=repeat(pi_vals, inner=length(size_vals)),
               mean_d=vec(mean_d),
               sd_d=vec(sd_d))
    return (summary=summary, effects=effects, beta_drift=beta,
            divergence=divergence, sizes=size_vals, pi_values=pi_vals)
end

"""
    qnca_resolution(X, Y; pi=1.0, n_grid=50, scope=nothing, monotone=:envelope)

Target-by-target account of one tolerance-indexed frontier. For each outcome
target `y` on the grid it reports the conditioning-set size `k`, the selected
type-1 rank, the raw order statistic, the fitted (monotone) requirement, whether
the fitted value was inherited from a lower target, the number and share of
attaining cases strictly below the fitted requirement, and the resistance
`rank - 1`: the number of added cases the raw ordinate absorbs without falling
below the smallest condition value of the original attaining cases.
"""
function qnca_resolution(X, Y; pi=1.0, n_grid=50, scope=nothing, monotone=:envelope)
    X = collect(float.(X)); Y = collect(float.(Y))
    validate_xy(X, Y)
    n_grid >= 2 || throw(ArgumentError("n_grid must be at least 2"))
    q = decimal_tail_probability(pi)
    mode = monotone_mode(monotone)
    x_min, x_max, y_min, y_max = validated_scope(X, Y, scope)
    y_grid = collect(range(y_min, y_max, length=n_grid))
    raw, k, h = raw_frontier(X, Y, q, y_grid)
    fitted = fit_frontier!(copy(raw), x_max, mode)
    m = length(y_grid)
    below = zeros(Int, m)
    share = fill(NaN, m)
    inherited = falses(m)
    for j in 1:m
        k[j] == 0 && continue
        below[j] = count(i -> Y[i] >= y_grid[j] && X[i] < fitted[j], eachindex(X))
        share[j] = below[j] / k[j]
        inherited[j] = fitted[j] > raw[j]
    end
    return (y=y_grid, k=k, rank=h, raw=raw, fitted=fitted, inherited=inherited,
            below=below, exception_share=share, resistance=max.(h .- 1, 0),
            pi=pi, monotone=mode)
end

function combinations_of(items::Vector{Int}, r::Int)
    out = Vector{Vector{Int}}()
    n = length(items)
    (r < 1 || r > n) && return out
    idx = collect(1:r)
    while true
        push!(out, items[idx])
        i = r
        while i >= 1 && idx[i] == n - r + i
            i -= 1
        end
        i == 0 && break
        idx[i] += 1
        for j in i+1:r
            idx[j] = idx[j-1] + 1
        end
    end
    return out
end

"""
    qnca_outliers(X, Y; pi_grid=[1.0, 0.95, 0.90], k=1, n_grid=50, scope=nothing,
                  monotone=:envelope, max_candidates=12, min_dif=0.01)

Deletion influence on `d_pi` across the tolerance grid, in the form of the NCA
outlier screen. With `k = 1` every case is deleted in turn. With `k > 1` all
combinations of `k` cases are deleted among the `max_candidates` cases lying
farthest below the most tolerant frontier at a target they attain, which keeps
cases that mask one another together. The scope stays fixed at the full-sample
rectangle. Rows are sorted by the largest absolute change over the grid, and a
row is flagged when a relative change reaches `min_dif` at some tolerance.
"""
function qnca_outliers(X, Y; pi_grid=[1.0, 0.95, 0.90], k=1, n_grid=50, scope=nothing,
                       monotone=:envelope, max_candidates=12, min_dif=0.01)
    X = collect(float.(X)); Y = collect(float.(Y))
    validate_xy(X, Y)
    n = length(X)
    n_grid >= 2 || throw(ArgumentError("n_grid must be at least 2"))
    (k >= 1 && k < n - 1) || throw(ArgumentError("k must be at least 1 and below n - 1"))
    pi_vals = collect(float.(pi_grid))
    isempty(pi_vals) && throw(ArgumentError("pi_grid must not be empty"))
    decimal_tail_probability.(pi_vals)
    mode = monotone_mode(monotone)
    bounds = validated_scope(X, Y, scope)
    d_full = [d_for_pi(X, Y, p, bounds, n_grid, mode) for p in pi_vals]

    if k == 1
        combos = [[i] for i in 1:n]
    else
        loose = pi_vals[argmin(pi_vals)]
        y_grid = collect(range(bounds[3], bounds[4], length=n_grid))
        phi = qnca_frontier_validated(X, Y, decimal_tail_probability(loose), y_grid;
                                      x_max=bounds[2], monotone=mode)
        gap = [phi[searchsortedlast(y_grid, Y[i])] - X[i] for i in 1:n]
        cand = [i for i in sortperm(gap; rev=true) if gap[i] >= 0.0]
        cand = cand[1:min(length(cand), max_candidates)]
        combos = combinations_of(cand, k)
    end

    d_without = Matrix{Float64}(undef, length(combos), length(pi_vals))
    @threads for r in eachindex(combos)
        keep = trues(n); keep[combos[r]] .= false
        for p in eachindex(pi_vals)
            d_without[r, p] = d_for_pi(X[keep], Y[keep], pi_vals[p], bounds, n_grid, mode)
        end
    end
    dif_abs = d_without .- reshape(d_full, 1, :)
    dif_rel = dif_abs ./ reshape(d_full, 1, :)
    dif_rel[:, d_full .== 0.0] .= NaN
    order = sortperm([maximum(abs.(dif_abs[r, :])) for r in eachindex(combos)]; rev=true)
    flagged = [any(x -> isfinite(x) && abs(x) >= min_dif, dif_rel[r, :]) for r in order]
    return (pi=pi_vals, d=d_full, cases=combos[order], d_without=d_without[order, :],
            dif_abs=dif_abs[order, :], dif_rel=dif_rel[order, :], flagged=flagged,
            k=k, monotone=mode)
end

"""
    qnca_sensitivity(X, Y; points, counts=0:5, pi_grid=[1.0, 0.95, 0.90],
                     n_grid=50, scope=nothing, monotone=:envelope)

Empirical breakdown curve. Adds `c` cases for each `c` in `counts` and
recomputes `d_pi` on the fixed scope of the original sample. `points` is one
`(x, y)` position, repeated `c` times, or a list of at least `maximum(counts)`
positions, of which the first `c` are added. Returns the effects and their
ratio to the effect without added cases.
"""
function qnca_sensitivity(X, Y; points, counts=0:5, pi_grid=[1.0, 0.95, 0.90],
                          n_grid=50, scope=nothing, monotone=:envelope)
    X = collect(float.(X)); Y = collect(float.(Y))
    validate_xy(X, Y)
    n_grid >= 2 || throw(ArgumentError("n_grid must be at least 2"))
    pi_vals = collect(float.(pi_grid))
    isempty(pi_vals) && throw(ArgumentError("pi_grid must not be empty"))
    decimal_tail_probability.(pi_vals)
    mode = monotone_mode(monotone)
    bounds = validated_scope(X, Y, scope)
    cs = collect(Int.(counts))
    (isempty(cs) || any(<(0), cs)) && throw(ArgumentError("counts must be non-negative"))
    pts = points isa Tuple ? [points] : collect(points)
    all(p -> length(p) == 2, pts) || throw(ArgumentError("each point must be an (x, y) pair"))
    for (px, py) in pts
        (bounds[1] <= px <= bounds[2] && bounds[3] <= py <= bounds[4]) ||
            throw(ArgumentError("added points must lie inside the scope"))
    end
    length(pts) == 1 || length(pts) >= maximum(cs) ||
        throw(ArgumentError("supply one point or at least maximum(counts) points"))
    d = Matrix{Float64}(undef, length(cs), length(pi_vals))
    for (r, c) in enumerate(cs)
        added = length(pts) == 1 ? fill(pts[1], c) : pts[1:c]
        Xc = vcat(X, Float64[float(p[1]) for p in added])
        Yc = vcat(Y, Float64[float(p[2]) for p in added])
        for p in eachindex(pi_vals)
            d[r, p] = d_for_pi(Xc, Yc, pi_vals[p], bounds, n_grid, mode)
        end
    end
    base = findfirst(==(0), cs)
    retention = base === nothing ? fill(NaN, size(d)) : d ./ reshape(d[base, :], 1, :)
    return (counts=cs, pi=pi_vals, d=d, retention=retention, monotone=mode)
end

end # module
