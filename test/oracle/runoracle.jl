# SPDX-License-Identifier: MPL-2.0
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# Agreement between the native methods and the Bioconductor references.
# Run: julia --project=test/oracle test/oracle/runoracle.jl
#
# Tolerance rule: both implementations are Monte-Carlo, so exact agreement is
# not expected. The yardstick is R against itself under a different seed. The
# native result must sit no further from R than SLACK × that seed-to-seed
# distance, plus a small absolute floor.

using CompositionalDA, DataFrames, RCall, Statistics, Test

include("datasets.jl")

const MC = 512
const SLACK = 2.0
const FLOOR = (effect = 0.05, pvalue = 0.01, wi_pvalue = 0.01)

"""Largest absolute difference of column `c` between two result tables."""
maxdiff(a, b, c) = maximum(abs.(a.table[!, c] .- b.table[!, c]))

"""
    agree(native, r1, r2) -> NamedTuple

Native-vs-R and R-vs-R(other seed) distances for effect, expected Welch p and
expected Wilcoxon p,
and whether native is within tolerance of R for each.
"""
function agree(native, r1, r2)
    out = map((:effect, :pvalue, :wi_pvalue)) do c
        dn = maxdiff(native, r1, c); dr = maxdiff(r2, r1, c)
        c => (native = dn, r_vs_r = dr, ok = dn <= SLACK * dr + FLOOR[c])
    end
    NamedTuple(out)
end

@testset "ALDEx2 oracle" begin
    datasets = Any[("mock", mock_dataset()), ("gut", fixture_dataset("gut")), ("soil", fixture_dataset("soil"))]
    for (name, (x, f, g, truth)) in datasets
        nat = aldex2(x, f, g; mc_samples = MC, seed = 1)
        r1 = aldex2_reference(x, f, g; mc_samples = MC, seed = 1)
        r2 = aldex2_reference(x, f, g; mc_samples = MC, seed = 2)
        a = agree(nat, r1, r2)
        println(name, ": ", a, "  R=", r1.provenance.parameters)
        @test a.effect.ok
        @test a.pvalue.ok
        @test a.wi_pvalue.ok
        @test cor(nat.table.effect, r1.table.effect) > 0.99
        if !isempty(truth)
            hit(t) = Set(t.table.feature[t.table.qvalue .< 0.05])
            println(name, ": native calls ", sort(collect(hit(nat))), "; R calls ", sort(collect(hit(r1))))
            @test hit(nat) == hit(r1)
        end
    end
end
