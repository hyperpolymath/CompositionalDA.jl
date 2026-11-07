# SPDX-License-Identifier: MPL-2.0
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>

"""
    bh_adjust(p) -> Vector{Float64}

Benjamini–Hochberg adjusted p-values, returned in the input order. This is the
step-up procedure, equivalent to R's `p.adjust(p, "BH")`. `NaN` entries stay
`NaN` and are excluded from the count of hypotheses `m`, as `p.adjust` excludes
`NA`. Values outside [0, 1] are refused.
"""
function bh_adjust(p::AbstractVector{<:Real})
    out = fill(NaN, length(p))
    idx = findall(!isnan, p)
    m = length(idx)
    m == 0 && return out
    # Validate before conversion: an invalid high-precision value just above
    # one can otherwise round to 1.0 and be silently accepted.
    all(v -> 0 <= v <= 1, p[idx]) || throw(ArgumentError("p-values must lie in [0, 1]"))
    q = Float64.(p[idx])
    ord = sortperm(q; rev = true)
    running = 1.0
    adj = similar(q)
    for (r, k) in enumerate(ord)
        rank = m - r + 1
        # Form the scale first: (p*m)/m can round below p even at rank m.
        # m/rank >= 1, so this multiplication preserves q >= p in Float64.
        running = min(running, (m / rank) * q[k])
        adj[k] = running
    end
    out[idx] = adj
    return out
end
