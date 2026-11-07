# SPDX-License-Identifier: MPL-2.0
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>

using CompositionalDA
using DataFrames
using Statistics
using Test
using Aqua

include("properties.jl")
include("mutants.jl")

@testset "CompositionalDA" begin
    @testset "Aqua" begin
        Aqua.test_all(CompositionalDA)
    end
    include("unit_core.jl")
    include("fuzz_core.jl")
    include("fuzz_statistics.jl")
    include("property_clr.jl")
    include("property_bh.jl")
    include("property_aldex2.jl")
    include("boundary_aldex2.jl")
    include("method_ancombc.jl")
    include("fuzz_ancombc.jl")
    include("posterior_aggregation.jl")
end
