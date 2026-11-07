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
using CompositionalDA: Provenance, DAResult, validate_counts, CompositionError,
                       validate_two_groups, validate_mc
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
    design = validate_two_groups(counts, features, groups)
    validate_mc(mc_samples, seed)
    seed <= typemax(Int32) || throw(ArgumentError("R reference seed must fit Int32"))
    levels = design.levels
    keep = vec(any(>(0), counts; dims = 1))
    count(keep) >= 2 || throw(CompositionError("ALDEx2 needs at least two nonzero features"))
    df = runner() do
        reads = permutedims(Float64.(counts[:, keep]))     # R integers are only 32-bit
        @rput reads
        conds = [g == levels[1] ? "0" : "1" for g in design.groups]; lv = ["0", "1"]
        @rput conds lv
        fnames = design.features[keep]; @rput fnames
        R"""
        set.seed($seed)
        reads <- as.data.frame(reads); rownames(reads) <- fnames
        colnames(reads) <- paste0("s", seq_len(ncol(reads)))
        # A character vector, not a factor: ALDEx2 1.42 calls round() on conds,
        # which errors for a factor. It builds its own factor, so the level
        # order is preserved by the explicit 0/1 encoding, even for original
        # numeric labels such as 2 and 10 whose lexical order differs.
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
                      contrast = "$(levels[2]) vs $(levels[1])", log_base = 2,
                      excluded_all_zero = design.features[.!keep])
    return DAResult(table, prov)
end

function CompositionalDA.ancombc_reference(counts::AbstractMatrix{<:Integer}, features::AbstractVector,
    groups::AbstractVector; prevalence::Real=0.1, pseudocount::Real=1.0, structural_zeros::Bool=false,
    negative_lower_bound::Bool=false, conserve::Bool=true, tolerance::Real=1e-7,
    max_iterations::Integer=1000, runner=f -> f())
    design = validate_two_groups(counts, features, groups)
    0 <= prevalence < 1 || throw(ArgumentError("prevalence must be in [0,1)"))
    isfinite(pseudocount) && pseudocount > 0 || throw(ArgumentError("pseudocount must be finite and positive"))
    isfinite(tolerance) && 0 < tolerance < 1 || throw(ArgumentError("tolerance must be in (0,1)"))
    max_iterations > 0 || throw(ArgumentError("max_iterations must be positive"))
    table, version = runner() do
        reads = permutedims(Float64.(counts))
        fnames = design.features
        conds = [g == design.levels[1] ? "reference" : "contrast" for g in design.groups]
        @rput reads fnames conds prevalence pseudocount structural_zeros negative_lower_bound conserve tolerance max_iterations
        R"""
        rownames(reads) <- fnames
        colnames(reads) <- paste0("s", seq_len(ncol(reads)))
        meta <- data.frame(condition=factor(conds, levels=c("reference", "contrast")),
                           row.names=colnames(reads))
        reference <- ANCOMBC::ancombc(data=reads, aggregate_data=reads,
            meta_data=meta, tax_level=NULL, formula="condition",
            p_adj_method="BH", pseudo=pseudocount, prv_cut=prevalence, lib_cut=0, group="condition",
            struc_zero=structural_zeros, neg_lb=negative_lower_bound, tol=tolerance,
            max_iter=max_iterations, conserve=conserve, global=FALSE, n_cl=1, verbose=FALSE)
        ref <- reference$res
        ids <- if ("taxon" %in% colnames(ref$lfc)) ref$lfc$taxon else rownames(ref$lfc)
        native_shape <- data.frame(feature=ids, effect=ref$lfc[["conditioncontrast"]],
            se=ref$se[["conditioncontrast"]], statistic=ref$W[["conditioncontrast"]],
            pvalue=ref$p_val[["conditioncontrast"]], qvalue=ref$q_val[["conditioncontrast"]])
        reference_version <- as.character(packageVersion("ANCOMBC"))
        """
        rcopy(DataFrame, R"native_shape"), rcopy(String, R"reference_version")
    end
    return DAResult(table, Provenance(:ancombc, :r_reference;
        r_package="ANCOMBC $version", prevalence, pseudocount, structural_zeros, negative_lower_bound,
        conserve, tolerance, max_iterations,
        contrast="$(design.levels[2]) vs $(design.levels[1])", log_base="e"))
end

end # module
