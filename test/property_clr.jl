# SPDX-License-Identifier: MPL-2.0
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>

@testset "P1 clr (mirrors proofs/agda CLR)" begin
    @testset "P1a clr sums to zero" begin
        @test prop_clr_sums_to_zero(clr)
        @test !prop_clr_sums_to_zero(mutant_clr_uncentred)
        @test !prop_clr_sums_to_zero(mutant_clr_first_ref)
    end
    @testset "P1b clr scale invariant" begin
        @test prop_clr_scale_invariant(clr)
        @test !prop_clr_scale_invariant(mutant_clr_uncentred)
    end
    @testset "P1c clr permutation equivariant" begin
        @test prop_clr_permutation_equivariant(clr)
        @test !prop_clr_permutation_equivariant(mutant_clr_position_weighted)
    end
end
