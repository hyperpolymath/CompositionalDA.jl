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

if haskey(ENV,"ALDEX2_REPAIR_DIR")
    repair_dir = ENV["ALDEX2_REPAIR_DIR"]
    @rput repair_dir
    R"""
    repair_environment <- new.env(parent=asNamespace("ALDEx2"))
    sys.source(file.path(repair_dir,"stats.fast.R"),envir=repair_environment)
    sys.source(file.path(repair_dir,"clr_ttest.r"),envir=repair_environment)
    for (symbol in c("t.fast","wilcox.fast","aldex.ttest"))
        assignInNamespace(symbol,repair_environment[[symbol]],ns="ALDEx2")
    for (test_file in c("test-directional-tails.R","test-aldex.ttest.R"))
        testthat::test_file(file.path(repair_dir,test_file),stop_on_failure=TRUE)
    """
    println("REFERENCE OVERRIDE: inclusive directional tails from ",repair_dir)
end

include("datasets.jl")
include("adapter_boundaries.jl")

const MC = 512
const SLACK = 2.0
const FLOOR = (effect = 0.05, pvalue = 0.01, qvalue = 0.01,
               wi_pvalue = 0.01, wi_qvalue = 0.01,
               diff_btw = 0.05, diff_win = 0.05, rab_all = 0.05)

"""Largest absolute difference of column `c` between two result tables."""
maxdiff(a, b, c) = maximum(abs.(a.table[!, c] .- b.table[!, c]))

"""
    agree(native, r1, r2) -> NamedTuple

Native-vs-R and R-vs-R(other seed) distances for effect, expected Welch p and
expected Wilcoxon p,
and whether native is within tolerance of R for each.
"""
function agree(native, r1, r2)
    out = map(keys(FLOOR)) do c
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
        for column in keys(FLOOR)
            @test a[column].ok
        end
        @test cor(nat.table.effect, r1.table.effect) > 0.99
        if !isempty(truth)
            hit(t) = Set(t.table.feature[t.table.qvalue .< 0.05])
            println(name, ": native calls ", sort(collect(hit(nat))), "; R calls ", sort(collect(hit(r1))))
            @test hit(nat) == hit(r1)
        end
    end
end

# The p-value and BH aggregations are deterministic conditional on the draws.
# Check identical R draws independently of both random-number generators and
# of the deliberately different effect-size pairing correlations.
@testset "ALDEx2 identical-instance oracle" begin
    for (name, (x, f, g, _)) in [("mock", mock_dataset()), ("gut", fixture_dataset("gut")), ("soil", fixture_dataset("soil"))]
        reads = permutedims(Float64.(x)); fnames = string.(f); conds = string.(g)
        @rput reads fnames conds
        R"""
        set.seed(73)
        rownames(reads) <- fnames
        colnames(reads) <- paste0("s", seq_len(ncol(reads)))
        clr_reference <- ALDEx2::aldex.clr(reads, conds, mc.samples=64, denom="all", verbose=FALSE)
        test_reference <- ALDEx2::aldex.ttest(clr_reference, verbose=FALSE)
        draws <- ALDEx2::getMonteCarloInstances(clr_reference)
        shared_draws <- array(unlist(draws), dim=c(nrow(reads), 64L, ncol(reads)))
        """
        inst = permutedims(rcopy(Array{Float64,3}, R"shared_draws"), (3, 1, 2))
        ref = rcopy(DataFrame, R"test_reference")
        levels = sort(unique(g))
        i1, i2 = (findall(==(level), g) for level in levels)
        native = CompositionalDA.aldex2_from_instances(inst, i1, i2;
            rng=StableRNG(73), first_is_i1=g[1] == levels[1])
        for (julia_col, r_col) in [(:pvalue, :we_ep), (:qvalue, :we_eBH),
                                  (:wi_pvalue, :wi_ep), (:wi_qvalue, :wi_eBH)]
            gap = maximum(abs.(native[!, julia_col] .- ref[!, r_col]))
            println(name, " shared draws ", julia_col, " maximum error = ", gap)
            @test gap <= 1e-12
        end
    end
end
