# SPDX-License-Identifier: MPL-2.0
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# Planted defects. A property test only counts as evidence if it rejects the
# mutant aimed at it; each property_*.jl file asserts exactly that.

"""Mutant: log without centring. Breaks P1a and P1b."""
mutant_clr_uncentred(x) = log.(x)

"""Mutant: centre on the first feature instead of the mean (an ALR in disguise). Breaks P1a only."""
mutant_clr_first_ref(x) = log.(x) .- log.(x[:, 1:1])

"""Mutant: centre with a fixed weight on the last column. Breaks P1c only."""
mutant_clr_position_weighted(x) = (l = log.(x); w = collect(1.0:size(x, 2))'; l .- sum(l .* w; dims = 2) ./ sum(w))

"""Mutant: Bonferroni-style p·m / rank with no running minimum. Breaks P5 monotonicity."""
function mutant_bh_no_cummin(p)
    out = fill(NaN, length(p)); idx = findall(!isnan, p); m = length(idx)
    ord = sortperm(p[idx]); q = similar(p[idx], Float64)
    for (r, k) in enumerate(ord); q[k] = min(1.0, p[idx][k] * m / r); end
    out[idx] = q; out
end

"""Mutant: instances are log-proportions without centring. Breaks P3."""
function mutant_instances_uncentred(rng, counts)
    inst = CompositionalDA.dirichlet_clr_instances(rng, counts; mc_samples = 8)
    inst .+ reshape(1:size(inst, 1), :, 1, 1)
end

"""Mutant: draws fresh effect-size pairings per feature, so the effect depends
on feature position. Breaks P6."""
function mutant_aldex2_per_feature_pairs(inst, i1, i2; rng)
    t = CompositionalDA.aldex2_from_instances(inst, i1, i2; rng)
    n1 = length(i1) * size(inst, 3); n2 = length(i2) * size(inst, 3)
    for j in 1:size(inst, 2)
        p = CompositionalDA.effect_pairs(rng, n1, n2)
        t.effect[j] = CompositionalDA.aldex2_effect(p, vec(inst[i1, j, :]), vec(inst[i2, j, :]))[1]
    end
    t
end
