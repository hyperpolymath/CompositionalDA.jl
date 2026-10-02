# SPDX-License-Identifier: MPL-2.0
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# Oracle datasets. `mock` is synthetic with a known answer. The public gut and
# soil tables are fetched by checksum at run time and their sources and licences
# are recorded in docs/DATASETS.md.

using Distributions: Multinomial, LogNormal
using StableRNGs

"""
    mock_dataset(; seed = 2026) -> (counts, features, groups, truth)

A two-group synthetic count table with planted differential features: 10
samples per group, 40 features with log-normal base abundances, features 1–3
up 8×, 4–5 down 4× in group "b", and multinomial sampling at depths 5k–20k.
`truth` lists the planted features. It is not a MetaManifold fixture (that
repository has no count-table fixtures).
"""
function mock_dataset(; seed::Integer = 2026)
    r = StableRNG(seed)
    d = 40; n = 10
    base = rand(r, LogNormal(0, 1.5), d)
    fold = ones(d); fold[1:3] .= 8; fold[4:5] .= 0.25
    rows = Vector{Vector{Int}}()
    for g in 1:2, _ in 1:n
        w = base .* (g == 2 ? fold : ones(d)) .* rand(r, LogNormal(0, 0.3), d)
        push!(rows, rand(r, Multinomial(rand(r, 5_000:20_000), w ./ sum(w))))
    end
    counts = permutedims(reduce(hcat, rows))
    return counts, ["mock$j" for j in 1:d], vcat(fill("a", n), fill("b", n)), ["mock$j" for j in 1:5]
end
