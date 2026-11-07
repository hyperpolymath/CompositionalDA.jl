<!-- SPDX-License-Identifier: CC-BY-SA-4.0 -->
<!-- SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk> -->

# MetaManifold issue 23: acceptance boundary

The requested contribution comprises native ALDEx2 and ANCOM-BC, independently
useful as a Julia package, followed by a narrow integration into the existing
Differential Abundance route and tab. Multinomial/Dirichlet-multinomial models
belong to issue 22 and CompositionalCounts.jl.

Release requires all of the following, independently:

1. Julia tests on the declared supported versions, including the WebUI's
   pinned Julia 1.12.5, and package-quality checks.
2. Pinned Agda/stdlib proof checking, axiom audit and negative controls.
3. ALDEx2 comparison of **every** returned statistic, plus identical-instance
   Welch/Wilcoxon p-value and expected-BH comparisons. Comparing only
   dimensionless effects and p-values cannot detect the wrong logarithm base.
4. ANCOM-BC conformance to the declared reference version and settings.
   A failure must remain visible; numerical thresholds cannot be widened
   after observing results to turn disagreement into success.
5. A registered Julia release, installable by name/UUID with registry pins.
   Joshua's issue comment explicitly rejects URL-pinned dependencies.
6. WebUI integration on its current main, a registry-resolved Manifest,
   no changes to renv.lock, and all required backend/frontend/server checks.
7. Signed, GitHub-verified commits, no merge commits, correct per-file SPDX,
   passing REUSE and a PR describing versions, validation and limitations.

A local development environment may use `Pkg.develop` to test an integration.
It is not a compliant released dependency or a substitute for registration.
Do not submit an upstream WebUI PR with a development path or invented
registry tree hash in its Manifest. Registration/maintainer acceptance are
external conditions, not facts that local tests can establish.

## Mathematical dependencies inspected

Mathesis's discovery atlas already points to this package's CLR/BH proofs.
Its README explicitly makes it an index rather than a proof authority; new
method-specific proofs belong here. CompositionalCounts.jl has a different
model and does not supply the missing ANCOM-BC estimator or calibration
proof. JanusKey concerns reversible file operations; Absolute Zero concerns
null effects/disclosures; Trope Checker verifies composition of declared IR
grades. None of those inspected interfaces establishes a DA sampling model
or FDR guarantee. Residual evidence types are useful for distinguishing
claims, assumptions and missing evidence, but adding them as a runtime
dependency would not prove the numerical implementation.

The useful mathematical additions are therefore local: PosteriorP's finite
minimum/sum inequality and BiasModel's cancellation/observability identities.
The exact rational and ring statements, the floating-point implementation,
the reference comparison and the statistical assumptions remain separate.
