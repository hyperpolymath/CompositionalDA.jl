# SPDX-License-Identifier: MPL-2.0
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>

using CompositionalDA, DataFrames, RCall, Test
include("datasets.jl")

# Fixed BEFORE examining the comparison results. These are numerical
# conformance thresholds, not confidence bounds or proofs of calibration.
const ANCOMBC_TOLERANCE = (effect=0.05, se=0.02, pvalue=0.02, qvalue=0.02)
@testset "ANCOM-BC Bioconductor conformance" begin
    for (name, (x, f, g, _)) in [("mock", mock_dataset()), ("gut", fixture_dataset("gut")), ("soil", fixture_dataset("soil"))]
        native = ancombc(x, f, g; structural_zeros=false)
        reference = ancombc_reference(x, f, g; structural_zeros=false)
        @test Set(native.table.feature) == Set(reference.table.feature)
        positions = indexin(native.table.feature, reference.table.feature)
        @test all(!isnothing, positions)
        moments = CompositionalDA.ancombc_moments(log.(Float64.(x) .+ 1),
            findall(==(sort(unique(g))[1]), g), findall(==(sort(unique(g))[2]), g))
        println(name, " bias diagnostic: native EM=", native.provenance.parameters[:bias],
            " WLS=", native.provenance.parameters[:bias_wls],
            " reference implied bias=", median(moments.difference .- reference.table[Int.(positions), :effect]))
        for column in keys(ANCOMBC_TOLERANCE)
            gap = maximum(abs.(native.table[!, column] .- reference.table[Int.(positions), column]))
            println(name, " ANCOM-BC ", column, " maximum error = ", gap)
            @test gap <= ANCOMBC_TOLERANCE[column]
        end
    end
end
