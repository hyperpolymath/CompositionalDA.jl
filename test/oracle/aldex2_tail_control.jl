# SPDX-License-Identifier: MPL-2.0
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
using RCall, Test

@testset "Pristine ALDEx2 discrete-tail defect control" begin
    @test rcopy(String,R"as.character(packageVersion('ALDEx2'))") == "1.42.0"
    R"""
    rank_data <- rbind(1:4,4:1)
    upper <- ALDEx2:::wilcox.fast(rank_data,c(1,1,2,2),FALSE)
    legacy <- pmin(1,2*upper,2*(1-upper))
    inclusive <- vapply(seq_len(2),function(i) wilcox.test(rank_data[i,1:2],
        rank_data[i,3:4],alternative="two.sided")$p.value,0)
    """
    @test rcopy(Vector{Float64},R"inclusive") ≈ [1/3,1/3]
    @test rcopy(Float64,R"legacy[1]") < 1e-12
    @test rcopy(Float64,R"legacy[2]") ≈ 1/3
    println("KNOWN-DEFECT CONTROL: pristine ALDEx2 1.42.0 drops the observed rank atom; legacy values = ",rcopy(R"legacy"))
end
