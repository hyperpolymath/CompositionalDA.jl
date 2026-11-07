# SPDX-License-Identifier: MPL-2.0
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# ALDEx2, written from docs/ALDEx2-SPEC.md and Fernandes et al. 2013 (PLoS ONE
# 8:e67019) and 2014 (Microbiome 2:15). It is not a port of the Bioconductor source.

"""
    dirichlet_clr_instances(rng, counts; mc_samples = 128, prior = 0.5) -> Array{Float64,3}

Monte-Carlo instances of the CLR. For each sample i and instance k, draw
proportions from Dirichlet(counts[i, :] .+ prior) by normalising Gamma(α, 1)
draws, then take the CLR of the draw. The result is indexed
`[sample, feature, instance]`.

Every instance is centred: each sample's row sums to zero (tested as P3 in
`test/property_aldex2.jl`, not proved). The instances are deliberately not
depth-invariant: a deeper sample has a tighter Dirichlet posterior, which is how
ALDEx2 carries sampling uncertainty into the tests.
"""
function dirichlet_clr_instances(rng::AbstractRNG, counts::AbstractMatrix{<:Integer};
                                 mc_samples::Integer = 128, prior::Real = 0.5)
    validate_counts(counts)
    mc_samples >= 2 || throw(ArgumentError("mc_samples must be at least 2"))
    isfinite(prior) && prior > 0 || throw(ArgumentError("prior must be finite and positive"))
    stored_prior = Float64(prior)
    isfinite(stored_prior) && stored_prior > 0 || throw(ArgumentError("prior must be representable as a positive finite Float64"))
    shapes = float_count_matrix(counts) .+ stored_prior
    all(isfinite, shapes) || throw(CompositionError("Dirichlet shape is not representable in Float64"))
    n, d = size(counts)
    out = Array{Float64,3}(undef, n, d, mc_samples)
    g = Vector{Float64}(undef, d)
    for i in 1:n, k in 1:mc_samples
        for j in 1:d
            g[j] = rand(rng, Gamma(shapes[i, j], 1.0))
        end
        # Flooring a zero draw would silently change the posterior. Refuse an
        # unrepresentable draw instead. ALDEx2 reports CLR values in log2 units.
        all(v -> isfinite(v) && v > 0, g) || throw(CompositionError("Dirichlet draw underflowed or overflowed; no result produced"))
        l = log2.(g)
        out[i, :, k] .= l .- mean(l)
    end
    return out
end

"""
    welch_tails(a, b) -> (upper, lower)

Both one-sided Welch t-test p-values, evaluated separately with
Welch–Satterthwaite degrees of freedom to preserve small probabilities.
"""
function welch_tails(a, b)
    variance = var(a) / length(a) + var(b) / length(b)
    isfinite(variance) && variance > 0 || throw(CompositionError("Welch test has zero or non-finite sampling variance"))
    try
        test = UnequalVarianceTTest(a, b)
        return pvalue(test; tail=:right), pvalue(test; tail=:left)
    catch error
        error isa DomainError && throw(CompositionError("Welch test is not estimable at the available numerical precision"))
        rethrow()
    end
end
welch_upper_p(a, b) = first(welch_tails(a, b))

"""
    mann_whitney_u(a, b) -> Float64

`#{(x, y) : x ∈ a, y ∈ b, x > y}`, with each tie counted as ½.
"""
function mann_whitney_u(a, b)
    u = 0.0
    for x in a, y in b
        u += x > y ? 1.0 : x == y ? 0.5 : 0.0
    end
    return u
end

"""
    mann_whitney_upper_tail(n1, n2) -> Vector{Float64}

Exact null upper tail of the Mann–Whitney statistic for group sizes `(n1, n2)`:
element `u + 1` is `P(U' ≥ u)` for `u = 0..n1·n2`.
"""
function mann_whitney_upper_tail(n1::Integer, n2::Integer)
    # prev[n+1] / cur[n+1] hold the point probabilities for sizes (m-1, n) / (m, n).
    prev = [[1.0] for _ in 0:n2]
    for m in 1:n1
        cur = Vector{Vector{Float64}}(undef, n2 + 1)
        cur[1] = [1.0]
        for n in 1:n2
            p = zeros(m * n + 1)
            for (u, v) in enumerate(prev[n + 1]); p[u + n] += v * m / (m + n); end
            for (u, v) in enumerate(cur[n]);      p[u]     += v * n / (m + n); end
            cur[n + 1] = p
        end
        prev = cur
    end
    return clamp.(reverse(cumsum(reverse(prev[n2 + 1]))), 0.0, 1.0)  # rounding can pass 1
end

"""
    wilcoxon_tails(a, b, tail) -> (upper, lower)

One-sided Wilcoxon rank-sum p-value for `a` above `b`. With `tail` from
[`mann_whitney_upper_tail`](@ref) it is exact, `P(U' ≥ U)`; with `tail ===
nothing` it is the normal approximation without continuity correction.
"""
function wilcoxon_tails(a, b, tail)
    u = mann_whitney_u(a, b)
    n1 = length(a); n2 = length(b)
    pooled = sort(vcat(a, b))
    all(isfinite, pooled) || throw(CompositionError("Wilcoxon observations must be finite"))
    correction = 0.0
    first = 1
    while first <= length(pooled)
        last = first
        while last < length(pooled) && pooled[last+1] == pooled[first]
            last += 1
        end
        size = Float64(last - first + 1)
        correction += size^3 - size
        first = last + 1
    end
    if tail !== nothing && correction == 0
        # Both tails include the observed atom. 1-upper is the STRICT lower
        # tail and gives direction-dependent, anti-conservative rank p-values.
        return tail[Int(u) + 1], tail[n1*n2 - Int(u) + 1]
    end
    n = n1 + n2
    variance = n1*n2 / 12 * (n+1 - correction / (n*(n-1)))
    isfinite(variance) && variance > 0 || throw(CompositionError("Wilcoxon test has zero or non-finite rank variance"))
    z = (u - n1*n2/2) / sqrt(variance)
    return ccdf(Normal(), z), ccdf(Normal(), -z)
end
wilcoxon_upper_p(a, b, tail) = first(wilcoxon_tails(a, b, tail))

"""
    effect_pairs(rng, n1, n2) -> NamedTuple

Index permutations for the ALDEx2 effect size over pooled groups of `n1` and
`n2` values, each truncated to `m = min(n1, n2)`: two per group for the
within-group spread (`w1a`, `w1b`, `w2a`, `w2b`), then one per group for the
between-group difference (`b1`, `b2`). They are drawn once and shared by every
feature, so relabelling the features permutes the effect sizes exactly (P6).
"""
function effect_pairs(rng::AbstractRNG, n1::Integer, n2::Integer)
    m = min(n1, n2)
    w1a = randperm(rng, n1)[1:m]; w1b = randperm(rng, n1)[1:m]
    w2a = randperm(rng, n2)[1:m]; w2b = randperm(rng, n2)[1:m]
    b1 = randperm(rng, n1)[1:m];  b2 = randperm(rng, n2)[1:m]
    return (; w1a, w1b, w2a, w2b, b1, b2)
end

"""
    aldex2_effect(pairs, x1, x2) -> (effect, diff_btw, diff_win)

ALDEx2 effect size for one feature. `x1` and `x2` are the instance-pooled CLR
values of the lower and higher sorted level, and `pairs` comes from
[`effect_pairs`](@ref). `btw = x2 − x1` over one pairing, `winmax` is the larger
of the two groups' within-group absolute differences, and the results are the
medians of `btw ./ winmax` (with `0/0` as 0), of `btw` and of `winmax`.
"""
function aldex2_effect(pairs, x1::AbstractVector, x2::AbstractVector)
    win = max.(abs.(x1[pairs.w1a] .- x1[pairs.w1b]), abs.(x2[pairs.w2a] .- x2[pairs.w2b]))
    btw = x2[pairs.b2] .- x1[pairs.b1]
    ratio = btw ./ win
    ratio[isnan.(ratio)] .= 0.0
    return median(ratio), median(btw), median(win)
end

"""
    aldex2_from_instances(inst, i1, i2; rng, first_is_i1 = true) -> DataFrame

The deterministic part of ALDEx2, applied to Monte-Carlo CLR instances
`inst[sample, feature, instance]`, with `i1`/`i2` the rows of the lower/higher
sorted level. Test group 1 is `i1` if `first_is_i1`, else `i2`. Each instance
gives one-sided Welch and Wilcoxon p-values for test group 1 above group 2; the
reported `pvalue`/`wi_pvalue` and `qvalue`/`wi_qvalue` are the consistent-sign
expected p and expected BH of `docs/ALDEx2-SPEC.md`. The effect columns are on
the sorted levels; `rng` is used only for their index permutations. The columns
are not named here; [`aldex2`](@ref) adds the feature names.
"""
function aldex2_from_instances(inst::AbstractArray{<:Real,3}, i1::AbstractVector{<:Integer},
                               i2::AbstractVector{<:Integer}; rng::AbstractRNG,
                               first_is_i1::Bool = true)
    _, d, mc = size(inst)
    t1, t2 = first_is_i1 ? (i1, i2) : (i2, i1)
    tail = length(t1) < 50 && length(t2) < 50 ? mann_whitney_upper_tail(length(t1), length(t2)) : nothing
    acc = Dict(k => zeros(d) for k in (:wg, :wl, :wgq, :wlq, :rg, :rl, :rgq, :rlq))
    pw = Vector{Float64}(undef, d); pr = Vector{Float64}(undef, d)
    lw = similar(pw); lr = similar(pr)
    for k in 1:mc
        for j in 1:d
            a = inst[t1, j, k]; b = inst[t2, j, k]
            pw[j], lw[j] = welch_tails(a, b)
            pr[j], lr[j] = wilcoxon_tails(a, b, tail)
            all(p -> isfinite(p) && 0 <= p <= 1, (pw[j], lw[j])) ||
                throw(CompositionError("Welch test is not estimable for feature column $j at Monte Carlo instance $k"))
            all(p -> isfinite(p) && 0 <= p <= 1, (pr[j], lr[j])) ||
                throw(CompositionError("Wilcoxon test is not estimable for feature column $j at Monte Carlo instance $k"))
        end
        for (greater, less, g, l, gq, lq) in ((pw, lw, :wg, :wl, :wgq, :wlq), (pr, lr, :rg, :rl, :rgq, :rlq))
            up = min.(1.0, 2 .* greater); dn = min.(1.0, 2 .* less)
            acc[g] .+= up; acc[l] .+= dn
            acc[gq] .+= bh_adjust(up); acc[lq] .+= bh_adjust(dn)
        end
    end
    expected(g, l) = min.(acc[g], acc[l]) ./ mc
    we_p = expected(:wg, :wl); we_q = expected(:wgq, :wlq)
    wi_p = expected(:rg, :rl); wi_q = expected(:rgq, :rlq)

    pairs = effect_pairs(rng, length(i1) * mc, length(i2) * mc)
    eff = zeros(d); btw = zeros(d); win = zeros(d); rab = zeros(d)
    for j in 1:d
        x1 = vec(inst[i1, j, :]); x2 = vec(inst[i2, j, :])
        eff[j], btw[j], win[j] = aldex2_effect(pairs, x1, x2)
        rab[j] = median(vec(inst[:, j, :]))
    end
    return DataFrame(effect = eff, pvalue = we_p, qvalue = we_q, wi_pvalue = wi_p,
                     wi_qvalue = wi_q, diff_btw = btw, diff_win = win, rab_all = rab)
end

"""
    aldex2(counts, features, groups; mc_samples = 128, seed = 1) -> DAResult

Two-group ALDEx2 on raw counts (samples × features). `groups` gives each
sample's condition and must have exactly two levels; the effect compares the
second level in sorted order against the first. The tests compare the first
sample's condition against the other. The result has the columns `feature`,
`effect`, `pvalue`/`qvalue` (consistent-sign expected Welch p and expected BH),
`wi_pvalue`/`wi_qvalue` (Wilcoxon), `diff_btw`, `diff_win` and `rab_all`.

Deterministic for a fixed `seed` (StableRNG). Numerical agreement is tested
against the Bioconductor reference. This is not a proof of calibration.
"""
function aldex2(counts::AbstractMatrix{<:Integer}, features::AbstractVector,
                groups::AbstractVector; mc_samples::Integer = 128, seed::Integer = 1)
    design = validate_two_groups(counts, features, groups)
    validate_mc(mc_samples, seed)
    (; levels, i1, i2) = design
    first_is_i1 = design.groups[1] == levels[1]
    keep = vec(any(>(0), counts; dims = 1))
    count(keep) >= 2 || throw(CompositionError("ALDEx2 needs at least two nonzero features"))

    rng = StableRNG(seed)
    inst = dirichlet_clr_instances(rng, counts[:, keep]; mc_samples)
    table = aldex2_from_instances(inst, i1, i2; rng, first_is_i1)
    insertcols!(table, 1, :feature => design.features[keep])
    prov = Provenance(:aldex2, :native; seed, mc_samples, test = "welch+wilcoxon",
                      contrast = "$(levels[2]) vs $(levels[1])", prior = 0.5,
                      log_base = 2, directional_tails = "inclusive; tie-adjusted normal for tied ranks",
                      excluded_all_zero = design.features[.!keep])
    return DAResult(table, prov)
end

"""
    aldex2_reference(counts, features, groups; mc_samples = 128, seed = 1, runner = f -> f())

The Bioconductor ALDEx2 reference, with the same input and output shape as
[`aldex2`](@ref). It is defined by the RCall package extension and is available
only after `using RCall` with ALDEx2 installed in R. `runner` wraps every R call,
so a host can serialise access to its R session. It is a test oracle and
transition aid; see `docs/R-SUNSET.md`.
"""
function aldex2_reference end
