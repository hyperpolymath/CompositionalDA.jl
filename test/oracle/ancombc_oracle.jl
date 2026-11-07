# SPDX-License-Identifier: MPL-2.0
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>

using CompositionalDA, DataFrames, RCall, Statistics, Test
include("datasets.jl")
include("adapter_boundaries.jl")

# Diagnostic comparison to an explicitly pinned upstream repair. The ordinary
# job leaves this unset and still tests the unmodified Bioconductor release.
if haskey(ENV, "ANCOMBC_ESTIMATOR_SOURCE")
    estimator_source = ENV["ANCOMBC_ESTIMATOR_SOURCE"]
    @rput estimator_source
    R"""
    estimator_environment <- new.env(parent=asNamespace("ANCOMBC"))
    sys.source(estimator_source, envir=estimator_environment)
    assignInNamespace(".bias_em", estimator_environment$.bias_em, ns="ANCOMBC")
    """
    println("REFERENCE OVERRIDE: bias estimator loaded from ", estimator_source)
end

# Fixed BEFORE examining the comparison results. These are numerical
# conformance thresholds, not confidence bounds or proofs of calibration.
const ANCOMBC_TOLERANCE = (effect=0.05, se=0.02, pvalue=0.02, qvalue=0.02)
const PRISTINE_DEFECT_CONTROL = get(ENV, "ANCOMBC_PRISTINE_DEFECT_CONTROL", "0") == "1"
if PRISTINE_DEFECT_CONTROL
    haskey(ENV, "ANCOMBC_ESTIMATOR_SOURCE") && error("the pristine control cannot load a repair")
    @test rcopy(String, R"as.character(packageVersion('ANCOMBC'))") == "2.12.0"
    # Exact equal-variance optimum: max(residual^2 - known_variance, 0).
    # This negative control establishes why the original reference cannot
    # qualify the native estimator. An unexpected repair must change the gate.
    R"""
    boundary_objective <- function(t) -dnorm(1, 0, sqrt(0.99+t), log=TRUE)
    boundary_result <- nloptr::neldermead(0, boundary_objective, lower=0)$par
    """
    @test rcopy(Float64, R"boundary_result") == 0
    @test rcopy(Bool, R"boundary_objective(0.01) < boundary_objective(boundary_result)")
    println("KNOWN-DEFECT CONTROL: unmodified ANCOMBC 2.12.0; three soil comparisons are expected failures, not conformance passes")
end
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
                c(iterations=iterNum, residual=epsilon,
                  bias=delta_new, left=l1_new, right=l2_new,
                  extra_left=kappa1_new, extra_right=kappa2_new,
                  weight_left=pi1_new, weight_null=pi0_new, weight_right=pi2_new,
                  log_likelihood=sum(log(
                    pi0_new*dnorm(beta, delta_new, sqrt(nu0)) +
                    pi1_new*dnorm(beta, delta_new+l1_new, sqrt(nu0+kappa1_new)) +
                    pi2_new*dnorm(beta, delta_new+l2_new, sqrt(nu0+kappa2_new)))))
        }))
    """
    reference = try
        ancombc_reference(x, f, g; structural_zeros=false,
            max_iterations=REFERENCE_MAX_ITERATIONS)
    finally
        R"untrace('.bias_em', where=asNamespace('ANCOMBC'))"
    end
    diagnostics = rcopy(Matrix{Float64}, R"do.call(rbind, oracle_bias_fits)")
    println("Reference fit (iterations, residual, bias, left, right, extra left/right, weights left/null/right, log likelihood): ", diagnostics)
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
            " log likelihood=", native.provenance.parameters[:log_likelihood],
            " reference implied bias=", median(moments.difference .- reference.table[Int.(positions), :effect]))
        for column in keys(ANCOMBC_TOLERANCE)
            gap = maximum(abs.(native.table[!, column] .- reference.table[Int.(positions), column]))
            println(name, " ANCOM-BC ", column, " maximum error = ", gap)
            if PRISTINE_DEFECT_CONTROL && name == "soil" && column in (:effect, :pvalue, :qvalue)
                @test_broken gap <= ANCOMBC_TOLERANCE[column]
            else
                @test gap <= ANCOMBC_TOLERANCE[column]
            end
        end
    end
end
