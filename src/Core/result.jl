# SPDX-License-Identifier: MPL-2.0
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>

"""
    Provenance(method, backend, seed, parameters)

What produced a result: the method name (e.g. `:aldex2`), the backend
(`:native`, or `:r_reference` from the optional RCall extension), the RNG seed
(`nothing` for deterministic methods), the method parameters, and the package
version.
"""
struct Provenance
    method::Symbol
    backend::Symbol
    seed::Union{Nothing,UInt64}
    parameters::Dict{Symbol,Any}
    package_version::VersionNumber
end

"""
    Provenance(method, backend; seed = nothing, parameters...)

Build a `Provenance` and stamp it with this package's version.
"""
Provenance(method::Symbol, backend::Symbol; seed = nothing, parameters...) =
    Provenance(method, backend, seed === nothing ? nothing : UInt64(seed),
               Dict{Symbol,Any}(parameters), pkgversion(@__MODULE__))

"""
    DAResult(table, provenance)

A differential-abundance result. `table` is a `DataFrame` with one row per
feature. It always has the columns `feature`, `effect`, `pvalue` and `qvalue`
(BH-adjusted), plus any method-specific columns.
"""
struct DAResult
    table::DataFrame
    provenance::Provenance
    function DAResult(table::DataFrame, provenance::Provenance)
        required = ("feature", "effect", "pvalue", "qvalue")
        missingcols = filter(c -> !(c in names(table)), collect(required))
        isempty(missingcols) || throw(ArgumentError("DAResult table lacks columns $(missingcols)"))
        new(table, provenance)
    end
end
