# SPDX-License-Identifier: MPL-2.0
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# The R reference backend. It is a transition aid and a test oracle; see
# docs/R-SUNSET.md for when it stops being maintained. Nothing here holds R
# state: every R call runs inside the caller's `runner`, so a host application
# that serialises its single R session (e.g. MetaManifold's `with_r_lock`)
# passes that function in.

module CompositionalDARCallExt

using CompositionalDA
using CompositionalDA: Provenance, DAResult, validate_counts, CompositionError
using DataFrames
using RCall

"""
    aldex2_reference(counts, features, groups; mc_samples = 128, seed = 1, runner = f -> f())

Run Bioconductor ALDEx2 (`aldex(..., test = "t", effect = TRUE, denom = "all")`)
on the same input as [`CompositionalDA.aldex2`](@ref) and return a `DAResult` with
the same column names. Requires R with the ALDEx2 package installed.
"""
function CompositionalDA.aldex2_reference(counts::AbstractMatrix{<:Integer}, features::AbstractVector,
                                          groups::AbstractVector; mc_samples::Integer = 128,
                                          seed::Integer = 1, runner = f -> f())
    validate_counts(counts)
    size(counts, 2) == length(features) || throw(CompositionError("feature names do not match columns"))
    levels = sort(unique(string.(groups)))
    length(levels) == 2 || throw(CompositionError("aldex2 needs exactly 2 groups"))
    df = runner() do
        reads = permutedims(Matrix{Int}(counts))           # ALDEx2 wants features × samples
        @rput reads
        conds = string.(groups); lv = levels
        @rput conds lv
        fnames = string.(features); @rput fnames
        R"""
        set.seed($seed)
        reads <- as.data.frame(reads); rownames(reads) <- fnames
        colnames(reads) <- paste0("s", seq_len(ncol(reads)))
        # A character vector, not a factor: ALDEx2 1.42 calls round() on conds,
        # which errors for a factor. It builds its own factor, so the level
        # order is the same sorted order as `lv`.
        stopifnot(identical(sort(unique(conds)), lv))
        res <- ALDEx2::aldex(reads, conds, mc.samples = $mc_samples,
                             test = "t", effect = TRUE, denom = "all", verbose = FALSE)
        res <- res[fnames, ]
        res$feature <- rownames(res)
        aldex_version <- as.character(packageVersion("ALDEx2"))
        """
        rcopy(DataFrame, R"res"), rcopy(String, R"aldex_version")
    end
    res, ver = df
    # rcopy turns R's dotted column names (we.ep, diff.btw, ...) into underscores.
    table = DataFrame(feature = res[!, "feature"], effect = res[!, "effect"], pvalue = res[!, "we_ep"],
                      qvalue = res[!, "we_eBH"], wi_pvalue = res[!, "wi_ep"], wi_qvalue = res[!, "wi_eBH"],
                      diff_btw = res[!, "diff_btw"], diff_win = res[!, "diff_win"], rab_all = res[!, "rab_all"])
    prov = Provenance(:aldex2, :r_reference; seed, mc_samples, r_package = "ALDEx2 $ver",
                      contrast = "$(levels[2]) vs $(levels[1])")
    return DAResult(table, prov)
end

end # module
