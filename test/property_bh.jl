# SPDX-License-Identifier: MPL-2.0
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>

@testset "P5 bh monotone" begin
    @test prop_bh_monotone(bh_adjust)
    @test !prop_bh_monotone(mutant_bh_no_cummin)
end
