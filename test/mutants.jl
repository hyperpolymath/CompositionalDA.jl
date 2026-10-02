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
