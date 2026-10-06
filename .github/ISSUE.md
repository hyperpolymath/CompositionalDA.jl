<!-- SPDX-License-Identifier: MPL-2.0 -->
<!-- Closable issue: description discrepancy (ANCOM-BC named available vs planned) -->

# Issue: Description distinguishes implemented ALDEx2 from planned ANCOM-BC

**Status:** Open (closable when done — see checklist below).

**What was wrong:** The package description (`README.md` line 7, `Project.toml`,
`src/CompositionalDA.jl` docstring) named ANCOM-BC as available alongside ALDEx2,
while the README status note and `docs/` said it lands in v0.3 (planned, not yet
in `src/`).

**Fix applied:** All descriptions updated to distinguish implemented ALDEx2
(`src/Methods/ALDEx2.jl`, native `aldex2`) from planned ANCOM-BC (v0.3,
not yet in `src/`).

**Checklist (closable only when all checked):**

- [x] `README.md` description corrected.
- [x] `README.md` status block clarified.
- [x] `Project.toml` added description field.
- [x] `src/CompositionalDA.jl` docstring updated.
- [x] `docs/R-SUNSET.md` criterion text updated.
- [x] ~~`.github/workflows/agda-cicd.yml` added (catalogue reference).~~ Removed 2026-10-06: the catalogue workflow it called (`standards/.github/workflows/echidna-verify.yml`) declares no `workflow_call` trigger, so every run died at startup with zero jobs (run 37231048694). The `proofs` job in `.github/workflows/ci.yml` (`proofs/check.sh`, Debian-pinned Agda 2.6.4.3) is the repo's proof gate and is green.
- [x] This issue file created with resolution criteria.

**No axiomatic claims made.** All statements verified against the current tree
(`git show HEAD:README.md`, `git show HEAD:Project.toml`, `git show HEAD:src/CompositionalDA.jl`).
