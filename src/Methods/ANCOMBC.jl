# SPDX-License-Identifier: MPL-2.0
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>

# Independent implementation of the two-group model in Lin & Peddada (2020),
# doi:10.1038/s41467-020-17041-7, equations 6, 14, 19--22, 31, 33 and 35.
# This is not ANCOM-BC2 and does not implement covariates or repeated measures.

function ancombc_moments(y, i1, i2)
    d = size(y, 2)
    means = zeros(2, d)
    variances = zeros(2, d)
    offsets = zeros(size(y, 1))
    for (g, rows) in enumerate((i1, i2))
        block = y[rows, :]
        means[g, :] = vec(mean(block; dims = 1))
        offsets[rows] = vec(mean(block; dims = 2)) .- mean(block)
        residual = block .- offsets[rows] .- means[g:g, :]
        variances[g, :] = vec(sum(abs2, residual; dims = 1)) ./ length(rows)^2
    end
    return (; difference = vec(means[2, :] .- means[1, :]),
            variance = vec(sum(variances; dims = 1)), offsets)
end

# Bounded one-dimensional minimisation. The variance objective need not be
# unimodal, so bracket every sampled local minimum, retaining the best result.
function variance_component(x, v, weights, centre)
    residual2 = (x .- centre).^2
    upper = maximum(residual2)
    upper == 0 && return 0.0
    objective(t) = sum(weights .* (log.(v .+ t) .+ residual2 ./ (v .+ t)))
    grid = vcat(0.0, exp.(range(log(minimum(v) * 1e-6), log(max(upper, minimum(v) * 1e-6)); length = 65)))
    values = objective.(grid)
    best = argmin(values)
    answer, score = grid[best], values[best]
    golden = (sqrt(5.0) - 1) / 2
    for j in 2:length(grid)-1
        values[j] <= values[j-1] && values[j] <= values[j+1] || continue
        lo, hi = grid[j-1], grid[j+1]
        a, b = hi - golden * (hi-lo), lo + golden * (hi-lo)
        fa, fb = objective(a), objective(b)
        for _ in 1:80
            if fa < fb
                hi, b, fb = b, a, fa
                a = hi - golden * (hi-lo); fa = objective(a)
            else
                lo, a, fa = a, b, fb
                b = lo + golden * (hi-lo); fb = objective(b)
            end
        end
        candidate = (lo + hi) / 2
        value = objective(candidate)
        if value < score
            answer, score = candidate, value
        end
    end
    return answer
end

function mixture_responsibilities(x, v, centres, extra, mixing)
    logweights = hcat([log(mixing[k]) .- 0.5 .* (log.(2π .* (v .+ extra[k])) .+
                       (x .- centres[k]).^2 ./ (v .+ extra[k])) for k in 1:3]...)
    maxima = maximum(logweights; dims = 2)
    weights = exp.(logweights .- maxima)
    total = sum(weights; dims = 2)
    return weights ./ total, sum(maxima .+ log.(total))
end

function ancombc_wls(x, v, centres, extra, mixing)
    # The reference WLS approximation assigns lower/upper tails by the fitted
    # component proportions. These are not maximum-posterior classifications.
    lower, upper = quantile(x, [mixing[1], 1 - mixing[3]])
    classes = [value < lower ? 1 : value >= upper ? 3 : 2 for value in x]
    weights = [1 / (v[i] + extra[classes[i]]) for i in eachindex(x)]
    shifts = centres .- centres[2]
    return (; wls = sum(weights .* (x .- shifts[classes])) / sum(weights),
            variance = 1 / sum(weights), classes)
end

function ancombc_bias(x, v; tolerance = 1e-7, max_iterations = 1000)
    length(x) >= 4 || throw(CompositionError("ANCOM-BC requires at least four estimable taxa"))
    all(z -> isfinite(z) && z > 0, v) || throw(CompositionError("ANCOM-BC has a zero or non-finite variance; no result produced"))
    centres = quantile(x, [0.1, 0.5, 0.9])
    extra = [var(x), 0.0, var(x)]
    mixing = [0.1, 0.8, 0.1]
    responsibilities, likelihood = mixture_responsibilities(x, v, centres, extra, mixing)
    for iteration in 1:max_iterations
        oldcentres = copy(centres)
        oldextra = copy(extra)
        oldmixing = copy(mixing)
        mixing = vec(mean(responsibilities; dims = 1))
        all(>(eps(Float64)), mixing) || throw(CompositionError("ANCOM-BC mixture component collapsed; no result produced"))
        precision = [responsibilities[:, k] ./ (v .+ extra[k]) for k in 1:3]
        masses = sum.(precision)
        centres = [sum(precision[k] .* x) / masses[k] for k in 1:3]
        # Weighted isotonic projection is the constrained maximiser in the
        # means for fixed variances: negative <= null <= positive.
        if centres[1] > centres[2]
            centres[1:2] .= sum(masses[1:2] .* centres[1:2]) / sum(masses[1:2])
        end
        if centres[2] > centres[3]
            centres[2:3] .= sum(masses[2:3] .* centres[2:3]) / sum(masses[2:3])
            if centres[1] > centres[2]
                centres .= sum(masses .* centres) / sum(masses)
            end
        end
        for k in (1, 3)
            candidate = variance_component(x, v, responsibilities[:, k], centres[k])
            objective(t) = sum(responsibilities[:, k] .* (log.(v .+ t) .+ (x .- centres[k]).^2 ./ (v .+ t)))
            objective(candidate) <= objective(extra[k]) && (extra[k] = candidate)
        end
        next_resp, next_likelihood = mixture_responsibilities(x, v, centres, extra, mixing)
        isfinite(next_likelihood) || throw(CompositionError("ANCOM-BC mixture is non-finite"))
        next_likelihood + 1e-9 * (1 + abs(likelihood)) >= likelihood ||
            throw(CompositionError("ANCOM-BC likelihood decreased; no result produced"))
        change = maximum(abs.(vcat(centres-oldcentres, extra-oldextra, mixing-oldmixing)))
        if abs(next_likelihood - likelihood) <= tolerance * (1 + abs(likelihood)) &&
                change <= tolerance * (1 + maximum(abs, vcat(centres, extra, mixing)))
            uncertainty = ancombc_wls(x, v, centres, extra, mixing)
            return (; bias = centres[2], wls = uncertainty.wls, variance = uncertainty.variance,
                    iterations = iteration, likelihood = next_likelihood,
                    mixing, centres, extra)
        end
        responsibilities, likelihood = next_resp, next_likelihood
    end
    throw(CompositionError("ANCOM-BC mixture did not converge in $max_iterations iterations; no result produced"))
end

"""
    ancombc(counts, features, groups; prevalence = 0.1, pseudocount = 1.0,
            structural_zeros = false, negative_lower_bound = false,
            conserve = true, tolerance = 1e-7, max_iterations = 1000)

Two-group ANCOM-BC (the original 2020 model), on samples × taxa raw counts.
The effect is a natural-log difference, sorted group 2 minus sorted group 1.
BH adjusts the Wald p-values of tested taxa only. Structural-zero decisions
are reported separately with `status` and undefined (`NaN`) numerical tests.
All-zero and below-prevalence taxa are omitted and named in provenance.

This estimator assumes enough non-differential taxa to identify sampling bias.
It does not make absolute abundance identifiable from counts alone. Calibration,
floating-point execution and equivalence to the R package are not formal proofs.
See `docs/ANCOMBC-SPEC.md` for the exact model and verification boundaries.
"""
function ancombc(counts::AbstractMatrix{<:Integer}, features::AbstractVector,
                 groups::AbstractVector; prevalence::Real = 0.1,
                 pseudocount::Real = 1.0, structural_zeros::Bool = false,
                 negative_lower_bound::Bool = false, conserve::Bool = true,
                 tolerance::Real = 1e-7, max_iterations::Integer = 1000)
    design = validate_two_groups(counts, features, groups)
    isfinite(prevalence) && 0 <= prevalence < 1 || throw(ArgumentError("prevalence must be in [0,1)"))
    isfinite(pseudocount) && pseudocount > 0 || throw(ArgumentError("pseudocount must be finite and positive"))
    isfinite(tolerance) && 0 < tolerance < 1 || throw(ArgumentError("tolerance must be in (0,1)"))
    pseudocount = Float64(pseudocount)
    tolerance = Float64(tolerance)
    isfinite(pseudocount) && pseudocount > 0 || throw(ArgumentError("pseudocount must be representable as a positive finite Float64"))
    0 < tolerance < 1 || throw(ArgumentError("tolerance must be representable as a Float64 strictly between 0 and 1"))
    max_iterations > 0 || throw(ArgumentError("max_iterations must be positive"))
    negative_lower_bound && !structural_zeros && throw(ArgumentError("negative_lower_bound requires structural_zeros"))
    present = counts .> 0
    keep = vec(any(present; dims = 1)) .& (vec(mean(present; dims = 1)) .>= prevalence)
    kept = findall(keep)
    status = fill("tested", length(kept))
    structural = falses(length(kept), 2)
    if structural_zeros
        for (g, rows) in enumerate((design.i1, design.i2)), (j, col) in enumerate(kept)
            fraction = mean(present[rows, col])
            structural[j, g] = fraction == 0 || (negative_lower_bound &&
                fraction - 1.959963984540054 * sqrt(fraction * (1-fraction) / length(rows)) <= 0)
        end
        for j in eachindex(kept)
            if all(structural[j, :])
                status[j] = "structural_zero_both"
            elseif any(structural[j, :])
                status[j] = "structural_zero"
            end
        end
    end
    tested = findall(==("tested"), status)
    length(tested) >= 4 || throw(CompositionError("fewer than four taxa remain for ANCOM-BC bias estimation"))
    y = log.(Float64.(counts[:, kept[tested]]) .+ Float64(pseudocount))
    all(isfinite, y) || throw(CompositionError("log counts are not representable in Float64"))
    moments = ancombc_moments(y, design.i1, design.i2)
    fit = ancombc_bias(moments.difference, moments.variance; tolerance, max_iterations)
    d = length(kept)
    effect, se, statistic, p = (fill(NaN, d) for _ in 1:4)
    # The reference API's conserve option changes the uncertainty estimate,
    # not the fitted coefficient. WLS supplies the bias-variance approximation.
    effect[tested] = moments.difference .- fit.bias
    se[tested] = sqrt.(moments.variance) .+ (conserve ? sqrt(fit.variance) : 0.0)
    statistic[tested] = effect[tested] ./ se[tested]
    p[tested] = 2 .* ccdf.(Normal(), abs.(statistic[tested]))
    table = DataFrame(; feature = design.features[kept], effect, pvalue = p,
        qvalue = bh_adjust(p), se, statistic, status,
        structural_zero_group1 = structural[:, 1], structural_zero_group2 = structural[:, 2])
    provenance = Provenance(:ancombc, :native; model = "Lin-Peddada-2020-two-group",
        contrast = "$(design.levels[2]) vs $(design.levels[1])", log_base = "e",
        prevalence, pseudocount, structural_zeros, negative_lower_bound, conserve,
        tolerance, max_iterations, excluded_features = design.features[.!keep],
        bias = fit.bias, bias_wls = fit.wls, bias_variance = fit.variance,
        iterations = fit.iterations, log_likelihood = fit.likelihood,
        mixture_weights = fit.mixing, bias_classification = "mixture-weight-quantiles", p_adjust = "BH")
    return DAResult(table, provenance)
end

"""
    ancombc_reference(counts, features, groups; runner = f -> f(), kwargs...)

Bioconductor ANCOMBC's original `ancombc` function, through the optional RCall
extension. A numerical reference, never a fallback for the native method.
"""
function ancombc_reference end
