<!--
SPDX-License-Identifier: CC-BY-SA-4.0
SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
-->

# Proof status

Machine-checked proofs for CompositionalDA.jl, in **Agda 2.6.4.3** with
**agda-stdlib 2.1** (Debian 13's `agda-bin` and `agda-stdlib` packages), under
`--safe --without-K`. Three words are kept apart throughout this file:

- **proved** — type-checked by Agda with zero postulates, zero holes and
  `--safe` on every module (the axiom audit fails the build otherwise). A
  theorem may take *named assumptions* as hypotheses or module parameters; the
  table says which.
- **tested** — exercised by a Julia property test in `test/`, which is also
  shown to reject a planted mutant in `test/mutants.jl`.
- **not done** — neither of the above. Listed, not implied.

> **The toolchain is the Debian 13 pin, not a vendored release.** `proofs/lib.sh`
> resolves `/usr/bin/agda` before anything on `PATH`, checks the version is
> 2.6.4.3 and the library is `standard-library-2.1`, and refuses any other.

## Gates

| Gate | Command | Result |
| --- | --- | --- |
| Type-check every proof, audit, then self-test | `proofs/check.sh` | PASS — exit 0, no warnings |
| Axiom audit on its own | `proofs/tests/axiom-audit.sh` | PASS — 7/7 modules reachable from `CompositionalDA.All`, no postulates, no FFI, no unsound flags, no holes, all `--safe` |
| Gate self-test (does the gate reject anything?) | `proofs/tests/gate-selftest.sh` | PASS — see "Last verified" below for the count |

The CI job `proofs` in `.github/workflows/ci.yml` runs `proofs/check.sh` in a
digest-pinned `debian:13-slim` container with `agda-bin=2.6.4.3-1+b2` and
`agda-stdlib=2.1-4` installed from apt. `check.sh` type-checks, audits and
then runs the self-test; both scripts resolve the prover and the library
through `proofs/lib.sh`, so they cannot disagree about what they checked
against.

The exact commands, exit codes and counts of the most recent run are in
"Last verified" at the end of this file.

## What is proved

Every module is reachable from `proofs/agda/CompositionalDA/All.agda`.

### P1 — the centred log-ratio transform

Module `CompositionalDA.CLR`, parametrised by a commutative ring `R`, a
commutative monoid `M` (the compositions; `×ₘ` is Aitchison's perturbation),
a `LogHom L : M → R` and a `Reciprocal ρ` for the part count `n`. With
`clr x i = log (x i) − ι · Σⱼ log (x j)`:

| Theorem | Statement | Named assumptions | Julia test mirrored |
| --- | --- | --- | --- |
| `clr-sum-zero` | `Σᵢ clr x i ≈ 0#` | `R`, `L`, `ρ` | `test/property_clr.jl` "P1a clr sums to zero" (`prop_clr_sums_to_zero`; mutants `mutant_clr_uncentred`, `mutant_clr_first_ref`) |
| `clr-scale-invariant` | `clr (scale k x) i ≈ clr x i` for every `k : Pos` | `R`, `L`, `ρ` | "P1b clr scale invariant" (`prop_clr_scale_invariant`; mutant `mutant_clr_uncentred`) |
| `clr-permute` | `clr (rearrange (π ⟨$⟩ʳ_) x) i ≈ clr x (π ⟨$⟩ʳ i)` for every permutation `π` | `R`, `L`, `ρ` | "P1c clr permutation equivariant" (`prop_clr_permutation_equivariant`; mutant `mutant_clr_position_weighted`) |
| `clr-perturb` | `clr (perturb x y) i ≈ clr x i + clr y i` | `R`, `L`, `ρ` | — (lemma for P1b; not separately tested) |
| `clr-constant` | `clr (replicate n k) i ≈ 0#` | `R`, `L`, `ρ` | — (lemma for P1b) |
| `centre-sum-zero` | `Σ (centre v) ≈ 0#` | `R`, `ρ` | — (lemma for P1a) |

Module `CompositionalDA.CLR.Instance` shows the assumptions are satisfiable,
so the theorems above are not vacuous, and supplies the negative control:

| Theorem | Statement | Named assumptions |
| --- | --- | --- |
| `additive-log` | `LogHom +-*-commutativeRing +-0-commutativeMonoid` — the monoid `(ℕ, +)`, the ring `ℚᵘ`, `log a = mkℚᵘ (+ a) 0` (`a/1`) | none |
| `reciprocal` | `Reciprocal +-*-commutativeRing (suc m)` for every `m`, `ι = mkℚᵘ (+ 1) m` (`1/(m+1)`) | none |
| `clr-not-shift-invariant` | `¬ (CLR₂.clr (CLR₂.perturb x₀ y₀) 0F ≃ CLR₂.clr x₀ 0F)` at `n = 2` | none — this is why `clr-scale-invariant` needs the perturbing composition constant |

### P5 — Benjamini–Hochberg step-up

Module `CompositionalDA.BenjaminiHochberg.StepUp` models the loop in
`src/Core/bh.jl` on exact rationals `ℚᵘ`:

```
scale m k p = p · m/(k+1)                           -- q[k] * m / rank, rank = k+1
stepUp m run (p ∷ ps) = (run ⊓ scale m k p) ∷ stepUp m (run ⊓ scale m k p) ps
adjust {m} = stepUp m 1ℚᵘ                           -- running = 1.0, m = length
```

Theorems on the **sorted** sequence (what the loop sees):

| Theorem | Statement | Named assumptions |
| --- | --- | --- |
| `adjust-≤1` | `lookup (adjust ps) i ≤ 1` | none |
| `adjust-descending` | `Linked _≥_ (adjust ps)` — the output never rises along the sorted order | none |
| `adjust-≥` | `lookup ps i ≤ lookup (adjust ps) i` — **q ≥ p** | input descending; every p in `[0, 1]` |
| `adjust-monotone` | `lookup ps i ≤ lookup ps j → lookup (adjust ps) i ≤ lookup (adjust ps) j` — **p_i ≤ p_j ⇒ q_i ≤ q_j**, at any two positions, ties included | input descending; every p `≥ 0` |
| `stepUp-tie` | positions holding the same p-value as the head receive the head's adjusted value | descending, `≥ 0` |
| `scale-≥` | `p ≤ scale m k p` when `k+1 ≤ m` and `0 ≤ p` | as stated |
| `scale-antitone` | `scale m (k+1) p ≤ scale m k p` when `0 ≤ p` | as stated |
| `stepUp-≤run`, `stepUp-descending`, `stepUp-≥`, `stepUp-tie-from`, `stepUp-monotone`, `desc-lookup`, `1≤ratio`, `ratio-antitone` | the generalised forms the four `adjust-*` theorems are instances of | as stated in the module |

Theorems at the **original** indices (what `test/property_bh.jl` reads), in
the sub-module `Unsorted p π` where `π : Permutation m m` stands for
`ord = sortperm(q; rev = true)` and `bh k = lookup (adjust sorted) (π ⟨$⟩ˡ k)`
is the `adj[k] = running` write-back:

| Theorem | Statement | Named assumptions | Julia test mirrored |
| --- | --- | --- | --- |
| `bh-≤1` | `bh k ≤ 1` | none | `test/property_bh.jl` "P5 bh monotone" (`prop_bh_monotone`: `q .<= 1`) |
| `bh-≥` | `lookup p k ≤ bh k` | `Linked _≥_ sorted`; `∀ k → 0 ≤ p k`; `∀ k → p k ≤ 1` | same (`q .>= p .- 1e-15`); also the known-answer vector in `test/unit_core.jl` (`bh_adjust([0.01, 0.04, 0.03, 0.005, 0.5, 0.04, NaN, 0.2])`) |
| `bh-monotone` | `lookup p i ≤ lookup p j → bh i ≤ bh j` | `Linked _≥_ sorted`; `∀ k → 0 ≤ p k` | same (`q[i] <= q[j] + 1e-15` whenever `p[i] <= p[j]`); mutant `mutant_bh_no_cummin` |
| `sorted-lookup` | `lookup sorted (π ⟨$⟩ˡ k) ≡ lookup p k` | none | — |

Module `CompositionalDA.BenjaminiHochberg` keeps the earlier per-rank
scaling as a single `mkℚᵘ`:

| Theorem | Statement | Named assumptions | Julia test mirrored |
| --- | --- | --- | --- |
| `bhScale-numerator`, `bhScale-denominator` | the numerator is `M·n`, the denominator `(j+1)(d+1)` (by `refl`) | none | — |
| `monotone-in-family-size` | `M ≤ M′ → bhScale M j n d ≤ bhScale M′ j n d` | `NonNegative n` | — (no Julia test varies the family size) |

### Support

`CompositionalDA.Prelude` (exact-rational lemmas: `common-denominator-+`,
`whole-is-one`, `zero-over`, `numerator-monotone-≤`, `≡⇒≃`, `+-cong-≃ʳ`) and
`CompositionalDA.LogHom` (the two assumption records) carry no claims of
their own.

## Named assumptions

Everything a theorem above takes as given, and how each is discharged — or
not. "Discharged by Julia" means a line of Julia code establishes it at run
time and that line is tested, not proved.

| Assumption | Where it appears | Discharged by |
| --- | --- | --- |
| `LogHom R M` — `log (a ×ₘ b) ≈ log a + log b` | module parameter of `CLR` | Proved satisfiable by `CLR.Instance.additive-log`. That the real logarithm on positive reals is one is a fact about ℝ, **not formalised**. |
| `Reciprocal R n` — `n · ι ≈ 1` | module parameter of `CLR` | Proved satisfiable by `CLR.Instance.reciprocal` for every `n ≥ 1`. `n = 0` has no reciprocal and no composition. |
| `Linked _≥_ sorted` — π sorts `p` descending | hypothesis of `Unsorted.bh-≥`, `Unsorted.bh-monotone` | Julia: `ord = sortperm(q; rev = true)` in `bh_adjust`. **Not proved**; `sortperm` is not re-implemented in Agda (`residue/benjamini-hochberg.residue` R-BH-3). |
| `0 ≤ p ≤ 1` for every p-value | hypotheses of `adjust-≥`, `adjust-monotone` and the `Unsorted` theorems | Julia: `bh_adjust` raises `ArgumentError` for any finite value outside `[0, 1]` before the loop runs; tested by `@test_throws ArgumentError bh_adjust([0.5, 1.5])` in `test/unit_core.jl`. **Not proved** (R-BH-4). The hypotheses are load-bearing: without `0 ≤ p` the tie argument fails; without `p ≤ 1`, `running = 1` undercuts the largest p-value. |
| `NonNegative n` | `monotone-in-family-size` | A p-value numerator is a count. Stated, not derived. |
| Rank bound `k + 1 ≤ m` | internal to `stepUp-≥` | **Proved** at the call site (`adjust-≥` passes `≤-refl`; the recursion passes `<⇒≤`). Not an assumption of any exported theorem. |

## Tested, not proved

Each item here is covered by a Julia test and by nothing in `proofs/`.

- **`sortperm` and the write-back.** That `ord` is a permutation sorting the
  non-NaN p-values descending, and that `adj[k] = running` writes sorted
  position `r` to original index `ord[r]` — i.e. that the Julia `ord` is the
  Agda `π` and `adj[k]` is `bh k`. Tested by `prop_bh_monotone` at the
  original indices.
- **NaN exclusion.** `ℚᵘ` has no NaN. The Agda model is the loop over the
  NaN-free vector with `m` its length; that NaN entries are excluded, that `m`
  counts only the rest, and that NaN-in gives NaN-out are Julia and tested:
  the generator `pvectors` in `test/properties.jl` plants a NaN in roughly one
  vector in five, `prop_bh_monotone` asserts
  `all(isnan, q[findall(isnan, p)])`, and `test/unit_core.jl` asserts
  `all(isnan, bh_adjust([NaN, NaN]))` and a known-answer vector containing a
  NaN.
- **The `[0, 1]` refusal.** The `ArgumentError` in `bh_adjust` is what makes
  the Agda hypotheses true at run time. The guard is tested
  (`@test_throws ArgumentError bh_adjust([0.5, 1.5])` in `test/unit_core.jl`);
  the proofs assume its postcondition.
- **Float64 ↔ `ℚᵘ`.** `q[k] * m / rank` in Float64 is a rounding of
  `scale m k p`; `min` on floats is `_⊓_` on rationals only up to that
  rounding. IEEE-754 is not modelled. The `1e-15` slack in `prop_bh_monotone`
  (`q .>= p .- 1e-15`, `q[i] <= q[j] + 1e-15`) is the test's allowance for
  this, not a proved bound. Likewise the real-valued `clr` in `src/` against
  the parametric theorems: `test/property_clr.jl` checks numerically with a
  tolerance.
- **Ties in floating point.** `adjust-monotone` covers ties exactly on
  rationals; whether two Float64 p-values that are equal stay equal through
  `* m / rank` is a floating-point question and is tested only. The generator
  does plant ties: `pvectors` rounds to 1–3 decimal digits (up to 40 entries)
  and, in the vectors where it does not plant a NaN, copies `p[1]` over a
  random position.
- **`monotone-in-family-size`** is proved but has **no Julia mirror**: no
  property test varies `m`.
- **R reference agreement** (ALDEx2/ANCOMBC through the optional RCall
  extension) is a test oracle only.

## Not done

- **FDR control.** That the q-values control the false discovery rate at the
  nominal level is a theorem about a probability measure over a sampling
  model. Nothing in `proofs/` bears on it and nothing in `test/` tests it
  against a ground truth; the README says the same. R-BH-7.
- **`bhScale` non-negativity in the expanded-denominator form** (R-BH-1):
  open for `bhScale`; the product form `scale` has it via `scale-≥`.
  `scale-is-bhScale`, which would connect the two modules, is not written.
- **A verified sort** replacing the `Linked _≥_ sorted` assumption (R-BH-3).
- **Real numbers.** No `LogHom` instance over ℝ, no `exp`, no continuity, no
  simplex geometry (R-CLR-1, R-CLR-5).
- **Floating-point models** of either method (R-CLR-4, R-BH-6).

The per-obligation records, with status and the shape of a closing proof,
are in `proofs/residue/clr.residue` and
`proofs/residue/benjamini-hochberg.residue`.

## The gate's own evidence

`proofs/tests/gate-selftest.sh` mutates a throwaway copy of the tree and
requires the gate to reject each mutation; a mutation whose `sed` pattern no
longer matches is a failure, not a pass. Each audit control must also fail with
the finding its row names, so one audit check cannot stand in for another; there
is one control per check in `axiom-audit.sh` (seven). The controls:

| Control | What is planted | Must be rejected by |
| --- | --- | --- |
| clr: sum-to-zero claimed as sum-to-one | `≈ 0#` → `≈ 1#` in `clr-sum-zero` | type-checker |
| clr: centring removed | `clr x = map log x` | type-checker |
| bh: family-size monotonicity reversed | `≤ℚ` operands swapped | type-checker |
| stepup: running minimum dropped | `stepUp` emits `scale m k p` without `⊓` (mirrors `mutant_bh_no_cummin`) | type-checker |
| stepup: q ≥ p reversed | conclusion of `adjust-≥` flipped | type-checker |
| stepup: monotonicity reversed | conclusion of `adjust-monotone` flipped | type-checker |
| stepup: scale written as rank/m | `mkℚᵘ (+ m) k` → `mkℚᵘ (+ suc k) m` | type-checker |
| audit: module dropped from the gate entry | `StepUp` import deleted from `All.agda` | axiom audit |
| audit: postulate injected | `postulate cheat : ∀ {A : Set} → A` | axiom audit |
| audit: `--safe` removed | from `CLR.agda` | axiom audit |
| audit: FFI pragma injected | `{-# FOREIGN GHC … #-}` appended to `Prelude.agda` | axiom audit |
| audit: unsound flag added | `--type-in-type` beside `--safe` in `Prelude.agda` | axiom audit |
| audit: trustMe mentioned | a `trustMe` reference appended to `Prelude.agda` | axiom audit |
| audit: hole left in a module | `hole = {! !}` appended to `Prelude.agda` | axiom audit |
| reject/FakeCLR, reject/MissingReciprocal, reject/Postulate | the standing negative controls | type-checker, with exit 42 **and** the error named on each file's `-- EXPECT:` line |
| pristine tree | nothing | must be **accepted** |

## Reproducing

```sh
proofs/check.sh                     # type-check, audit, then the gate self-test
proofs/tests/gate-selftest.sh       # prove the gate rejects planted defects
proofs/tests/axiom-audit.sh         # the audit alone
```

`check.sh` needs Debian 13's `agda-bin` 2.6.4.3 and `agda-stdlib` 2.1
(`apt install agda agda-stdlib`); `lib.sh` fails, naming the version it saw,
if anything else is found first. Nothing is downloaded and nothing is written
into the toolchain's data directory. **An absent prover is a failure, never a
skip** — in CI and locally alike.

## Last verified

2026-10-02, after rebasing onto PR #2 (`320513d`, which added `StepUp` and
the 18-control self-test), one run of the CI interface with Debian 13's
packages (`agda-bin 2.6.4.3-1+b2`, `agda-stdlib 2.1-4`), interface files
deleted first so every module was checked from source:

```
proofs/check.sh                  → exit 0   (check: Agda version 2.6.4.3, standard-library-2.1;
                                             axiom-audit: 7 module(s) reachable from CompositionalDA.All, clean;
                                             CompositionalDA/All.agda type-checks;
                                             gate-selftest: 18/18 controls behaved correctly)
proofs/tests/gate-selftest.sh    → exit 0   (gate-selftest: 18/18 controls behaved correctly)
proofs/tests/axiom-audit.sh      → exit 0   (auditing 7 module(s); axiom-audit: clean)
```

Wall-clock 14:58:02Z → 15:00:00Z. The only source change the 2.1 library
needed was `_≡?_` → `_≟_` on ℤ in `Prelude.agda`; `StepUp.agda` needed
none. The 18 controls are the 14 mutations, the 3 `reject/` files and the
pristine tree listed above; each mutation is confirmed to have changed its
file before the gate runs, and each `reject/` file is refused with exit 42
and the error its `-- EXPECT:` line names. The same 18/18 was obtained on
Agda 2.7.0.1 with the stdlib development SHA at PR #2 (14:18Z–14:22Z).

The Julia side, same run:

```
julia --project -e 'using Pkg; Pkg.test()'   → exit 0   (Test Summary: CompositionalDA | Pass 43 | Total 43)
```

That run includes `test/property_clr.jl` (P1a/P1b/P1c) and
`test/property_bh.jl` (P5), and the mutant checks in `test/mutants.jl`.
