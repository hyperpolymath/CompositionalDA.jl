# SPDX-License-Identifier: MPL-2.0
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# ALDEx2, written from Fernandes et al. 2013 (PLoS ONE 8:e67019) and 2014
# (Microbiome 2:15). It is not a port of the Bioconductor source.

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
    mc_samples > 0 || throw(ArgumentError("mc_samples must be positive"))
    prior > 0 || throw(ArgumentError("prior must be positive"))
    n, d = size(counts)
    out = Array{Float64,3}(undef, n, d, mc_samples)
    g = Vector{Float64}(undef, d)
    for i in 1:n, k in 1:mc_samples
        for j in 1:d
            g[j] = rand(rng, Gamma(counts[i, j] + prior, 1.0))
        end
        # Gamma draws with tiny shape can underflow to 0; floor them at the
        # smallest positive normal so the log stays finite.
        g .= max.(g, floatmin(Float64))
        l = log.(g)
        out[i, :, k] .= l .- mean(l)
    end
    return out
end

"""
    welch_p(a, b) -> Float64

Two-sided Welch two-sample t-test p-value.
"""
welch_p(a, b) = pvalue(UnequalVarianceTTest(a, b))

"""
    wilcoxon_p(a, b) -> Float64

Two-sided Wilcoxon rank-sum (Mann–Whitney) p-value. It is exact when both groups
have fewer than 50 observations and there are no ties, and uses the normal
approximation otherwise. This is the same rule as R's `wilcox.test` default.
"""
function wilcoxon_p(a, b)
    exact = length(a) < 50 && length(b) < 50 && allunique(vcat(a, b))
    exact ? pvalue(ExactMannWhitneyUTest(a, b)) : pvalue(ApproximateMannWhitneyUTest(a, b))
end

"""
    effect_pairs(rng, n1, n2, m) -> NamedTuple

Random index vectors for the ALDEx2 effect size: `m` draws of a group-1 element,
a group-2 element, and a second element of each group for the within-group
spread. They are drawn once and shared by every feature, so relabelling the
features permutes the effect sizes exactly (P6).
"""
effect_pairs(rng::AbstractRNG, n1::Integer, n2::Integer, m::Integer) =
    (a = rand(rng, 1:n1, m), b = rand(rng, 1:n2, m),
     a2 = rand(rng, 1:n1, m), b2 = rand(rng, 1:n2, m))

"""
    aldex2_effect(pairs, x1, x2) -> (effect, diff_btw, diff_win)

ALDEx2 effect size for one feature. `x1` and `x2` are the instance-pooled CLR
values of the two groups and `pairs` comes from [`effect_pairs`](@ref). The
between-group difference is taken over random cross-group pairs. The
within-group dispersion is the larger of the two groups' absolute differences
over random within-group pairs. The effect is the median of their ratio;
`diff_btw` and `diff_win` are the medians of each part.
"""
function aldex2_effect(pairs, x1::AbstractVector, x2::AbstractVector)
    a = x1[pairs.a]; b = x2[pairs.b]
    btw = b .- a
    win = max.(abs.(a .- x1[pairs.a2]), abs.(b .- x2[pairs.b2]))
    ratio = btw ./ max.(win, eps())
    return median(ratio), median(btw), median(win)
end

"""
    aldex2_from_instances(inst, i1, i2; rng) -> DataFrame

The deterministic part of ALDEx2, applied to Monte-Carlo CLR instances
`inst[sample, feature, instance]` with group-1 rows `i1` and group-2 rows `i2`.
For every instance it runs a Welch t-test and a Wilcoxon test per feature and
BH-adjusts across features. The reported values are the means over instances.
`rng` is used only for the effect-size pairing indices. The columns are not
named here; [`aldex2`](@ref) adds the feature names.
"""
function aldex2_from_instances(inst::AbstractArray{<:Real,3}, i1::AbstractVector{<:Integer},
                               i2::AbstractVector{<:Integer}; rng::AbstractRNG)
    _, d, mc = size(inst)
    we_p = zeros(d); we_q = zeros(d); wi_p = zeros(d); wi_q = zeros(d)
    pw = Vector{Float64}(undef, d); pr = Vector{Float64}(undef, d)
    for k in 1:mc
        for j in 1:d
            a = inst[i1, j, k]; b = inst[i2, j, k]
            pw[j] = welch_p(a, b)
            pr[j] = wilcoxon_p(a, b)
        end
        we_p .+= pw; we_q .+= bh_adjust(pw)
        wi_p .+= pr; wi_q .+= bh_adjust(pr)
    end
    we_p ./= mc; we_q ./= mc; wi_p ./= mc; wi_q ./= mc

    n1 = length(i1) * mc; n2 = length(i2) * mc
    pairs = effect_pairs(rng, n1, n2, max(n1, n2))
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
sample's condition and must have exactly two levels; the second level in sorted
order is compared against the first. The result has the columns `feature`,
`effect`, `pvalue`/`qvalue` (expected Welch p and expected BH),
`wi_pvalue`/`wi_qvalue` (Wilcoxon), `diff_btw`, `diff_win` and `rab_all`.

Deterministic for a fixed `seed` (StableRNG). Calibration is tested against the
Bioconductor reference, not proved.
"""
function aldex2(counts::AbstractMatrix{<:Integer}, features::AbstractVector,
                groups::AbstractVector; mc_samples::Integer = 128, seed::Integer = 1)
    validate_counts(counts)
    n, d = size(counts)
    length(features) == d || throw(CompositionError("$(length(features)) feature names for $d columns"))
    length(groups) == n || throw(CompositionError("$(length(groups)) group labels for $n samples"))
    levels = sort(unique(groups))
    length(levels) == 2 || throw(CompositionError("aldex2 needs exactly 2 groups, got $(length(levels))"))
    i1 = findall(==(levels[1]), groups)
    i2 = findall(==(levels[2]), groups)
    (length(i1) >= 2 && length(i2) >= 2) || throw(CompositionError("each group needs at least 2 samples"))

    rng = StableRNG(seed)
    inst = dirichlet_clr_instances(rng, counts; mc_samples)
    table = aldex2_from_instances(inst, i1, i2; rng)
    insertcols!(table, 1, :feature => collect(features))
    prov = Provenance(:aldex2, :native; seed, mc_samples, test = "welch+wilcoxon",
                      contrast = "$(levels[2]) vs $(levels[1])", prior = 0.5)
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
