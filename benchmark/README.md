<!-- SPDX-License-Identifier: CC-BY-SA-4.0 -->
<!-- SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) -->

# Exploratory launch measurements

From the repository root, prepare the separate environment:

```sh
julia --project=benchmark -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=benchmark benchmark/calibration_pilot.jl /path/to/pilot-output
julia --project=benchmark benchmark/benchmarks.jl /path/to/benchmark-output
```

Read `docs/CALIBRATION-PILOT-PROTOCOL.adoc` before interpreting the output.
Seeds, scenarios and comparisons were fixed before results. Refusals are
separate from successful fits. Neither script certifies calibration, speed,
<!-- coderabbit-ai-generated:v1 -->

<!-- coderabbit-ai-generated-section:v1:2daa59cb87a49803898dd374077604af362283cb0f06dd91a39c6483729ff4fc -->

stability or package readiness. Keep generated measurements outside source.
Each output directory includes the active environment's `Project.toml` and
resolved `Manifest.toml`, including transitive dependencies such as StatsFuns.
Retain both with the results; instantiate the active project before running.
Timing benchmarks refuse a known dirty source tree (including untracked files).
When Git provenance is unavailable, scripts warn and record `source_head` and
`source_tree` as `unknown`; those runs need a separate source archive to reproduce.

<!-- /coderabbit-ai-generated:v1 -->
Benchmark timings exclude compilation/warmup and retain trial counts; ten
trials provide an exploratory baseline, not precise tail-latency estimates.

The optional CI investigation exports the small-unbalanced complete-null
inputs and native results, then runs pinned ANCOMBC 2.12.0 with the explicitly
identified variance-estimator repair. It records R convergence rather than
treating an iteration-cap result as a valid fit. A successful investigation
job means the experiment completed, not that the adverse calibration passed.
<!-- coderabbit-ai-generated:v1 -->

<!-- coderabbit-ai-generated-section:v1:293a501349b835ccb7aa6a0ef06715bbeee5f83376f52e981bd78fa617317411 -->

The workflow requests 30 days of artifact retention; repository policy may
shorten that period. Reviewers should check the actual expiry before archiving
the raw artifact with its SHA-256 for a release decision.

<!-- /coderabbit-ai-generated:v1 -->
