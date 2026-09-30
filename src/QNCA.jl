"""
    QNCA

Quantile Necessary Condition Analysis: a tolerance-parameterised generalisation
of Dul's Necessary Condition Analysis. The deterministic NCA CE-FDH ceiling is
the `pi == 1` limit of the QNCA frontier family. Lowering `pi` discards the most
efficient cases and raises the necessity floor, trading crispness for robustness.

The implementation uses the same pool-adjacent-violators isotonic fit,
inverse-empirical-CDF quantile, and carry-forward of the frontier to the scope
ceiling as the R and Python siblings. Configured decimal tolerances are converted
to exact rational ranks before frontier evaluation.
"""
module QNCA

using Random, Statistics, Printf, Base.Threads

export generate_reverse_L, isotonic_increasing, qnca_frontier, qnca_d,
       nca_ce_fdh_d, qnca, qnca_rank, quantile_type1_pi,
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
    qnca_frontier(X, Y, pi, y_grid; x_max=nothing)

Quantile necessity frontier phi_pi(y), isotonic-fit and carried forward over
trailing (high-outcome) NaN levels. Under a fixed scope, pass `x_max` so the
whole high-outcome band is counted as empty, matching CE-FDH.
"""
function qnca_frontier(X, Y, pi, y_grid; x_max=nothing)
    Xf = collect(float.(X)); Yf = collect(float.(Y)); grid = collect(float.(y_grid))
    validate_xy(Xf, Yf)
    isempty(grid) && throw(ArgumentError("y_grid must not be empty"))
    all(isfinite, grid) || throw(ArgumentError("y_grid must contain only finite values"))
    issorted(grid) || throw(ArgumentError("y_grid must be non-decreasing"))
    q = decimal_tail_probability(pi)
    return qnca_frontier_validated(Xf, Yf, q, grid; x_max=x_max)
end

function qnca_frontier_validated(X, Y, q::Rational{BigInt}, y_grid; x_max=nothing)
    m = length(y_grid)
    phi = fill(NaN, m)
    @inbounds for j in 1:m
        idx = findall(>=(y_grid[j]), Y)
        if !isempty(idx)
            x = sort(collect(@view(X[idx])))
            rank = clamp(ceil(Int, length(x) * q), 1, length(x))
            phi[j] = x[rank]
        end
    end
    fin = findall(!isnan, phi)
    if length(fin) >= 2
        phi[fin] = isotonic_increasing(phi[fin])
    end
    if !isempty(fin)
        last = fin[end]
        if last < m
            phi[last+1:m] .= x_max === nothing ? phi[last] : x_max
        end
    end
    return phi
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

function d_for_pi(X, Y, pi, scope, n_grid)
    return qnca(X, Y; pi=pi, n_grid=n_grid, B=0, scope=scope).d_pi
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
              permutations=nothing)
    X = collect(float.(X)); Y = collect(float.(Y))
    validate_xy(X, Y)
    n_grid >= 2 || throw(ArgumentError("n_grid must be at least 2"))
    B >= 0 || throw(ArgumentError("B must be non-negative"))
    q = decimal_tail_probability(pi)
    x_min, x_max, y_min, y_max = validated_scope(X, Y, scope)
    scope_area = (x_max - x_min) * (y_max - y_min)
    y_grid = collect(range(y_min, y_max, length=n_grid))
    phi_obs = qnca_frontier_validated(X, Y, q, y_grid; x_max=x_max)
    d_obs = qnca_d(phi_obs, y_grid, x_min, scope_area)

    p_pi = NaN
    if permutations !== nothing
        perm = validate_permutations(permutations, length(Y))
        hits = zeros(UInt8, size(perm, 1))
        @threads for b in 1:size(perm, 1)
            Yp = Y[vec(perm[b, :])]
            d_p = qnca_d(qnca_frontier_validated(X, Yp, q, y_grid; x_max=x_max),
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
            d_p = qnca_d(qnca_frontier_validated(X, Yp, q, y_grid; x_max=x_max),
                          y_grid, x_min, scope_area)
            hits[b] = d_p >= d_obs ? 1 : 0
        end
        p_pi = (1 + sum(hits)) / (B + 1)
    end
    return (pi=pi, d_pi=d_obs, p_pi=p_pi, y_grid=y_grid, x_required=phi_obs, scope=scope_area)
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
                           uniform_draws=nothing)
    X = collect(float.(X)); Y = collect(float.(Y))
    validate_xy(X, Y)
    n_grid >= 2 || throw(ArgumentError("n_grid must be at least 2"))
    pi_vals = collect(float.(pi_grid))
    isempty(pi_vals) && throw(ArgumentError("pi_grid must not be empty"))
    decimal_tail_probability.(pi_vals)
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

    d_obs = [d_for_pi(X, Y, pi, fixed_scope, n_grid) for pi in pi_vals]
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
            d_null[j, k] = d_for_pi(Xj, Yj, pi_vals[k], fixed_scope, n_grid)
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
                           n_grid=50, scope=nothing, seed=nothing, subsamples=nothing)
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
                effects[k, rr, p] = d_for_pi(X[idx], Y[idx], pi_vals[p], fixed_scope, n_grid)
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

end # module
