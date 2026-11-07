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
stability or package readiness. Keep generated measurements outside source.
Benchmark timings exclude compilation/warmup and retain trial counts; ten
trials provide an exploratory baseline, not precise tail-latency estimates.

The optional CI investigation exports the small-unbalanced complete-null
inputs and native results, then runs pinned ANCOMBC 2.12.0 with the explicitly
identified variance-estimator repair. It records R convergence rather than
treating an iteration-cap result as a valid fit. A successful investigation
job means the experiment completed, not that the adverse calibration passed.
The raw artifact expires after 30 days; archive it with its SHA-256 before
using it in a release decision.
