# SPDX-License-Identifier: MPL-2.0
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>

using CompositionalDA, DataFrames, RCall, Statistics, Test
include("datasets.jl")

# Fixed BEFORE examining the comparison results. These are numerical
# conformance thresholds, not confidence bounds or proofs of calibration.
const ANCOMBC_TOLERANCE = (effect=0.05, se=0.02, pvalue=0.02, qvalue=0.02)
# ANCOMBC 2.12.0 silently returns its current parameters at max_iter. Its
# original 1000-iteration reference run disagreed on soil. Verify convergence
# explicitly; increasing the work budget is not relaxing numerical tolerances.
const REFERENCE_MAX_ITERATIONS = 20000
function checked_reference(x, f, g)
    R"""
    loadNamespace("ANCOMBC")
    oracle_bias_fits <- list()
    trace(".bias_em", where=asNamespace("ANCOMBC"), print=FALSE,
        exit=quote({
            .GlobalEnv$oracle_bias_fits[[length(.GlobalEnv$oracle_bias_fits)+1L]] <-
                c(iterations=iterNum, residual=epsilon)
        }))
    """
    reference = try
        ancombc_reference(x, f, g; structural_zeros=false,
            max_iterations=REFERENCE_MAX_ITERATIONS)
    finally
        R"untrace('.bias_em', where=asNamespace('ANCOMBC'))"
    end
    diagnostics = rcopy(Matrix{Float64}, R"do.call(rbind, oracle_bias_fits)")
    println("Reference convergence (iterations, parameter residual): ", diagnostics)
    @test size(diagnostics, 1) > 0
    @test all(diagnostics[:, 1] .< REFERENCE_MAX_ITERATIONS)
    @test all(diagnostics[:, 2] .<= 1e-7)
    return reference
end
@testset "ANCOM-BC Bioconductor conformance" begin
    for (name, (x, f, g, _)) in [("mock", mock_dataset()), ("gut", fixture_dataset("gut")), ("soil", fixture_dataset("soil"))]
        native = ancombc(x, f, g; structural_zeros=false)
        reference = checked_reference(x, f, g)
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
