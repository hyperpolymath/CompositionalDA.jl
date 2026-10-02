<!-- SPDX-License-Identifier: CC-BY-SA-4.0 -->
<!-- SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) -->

# Proof status

What is **proved** (machine-checked in Agda), what is only **tested** (Julia
property and unit tests), and what is still **open**. A row says only what its
column names; nothing else is implied.

Gate: Agda 2.7.0.1, `--safe --without-K`, no postulates, agda-stdlib pinned by SHA
(`proofs/bootstrap.sh`). `proofs/tests/axiom-audit.sh` checks every module is
`--safe`, has no postulate, FFI, unsound flag, trusted primitive or hole, and is
reachable from `CompositionalDA/All.agda`. `proofs/tests/gate-selftest.sh` shows
the gate rejects the three controls in `proofs/agda/reject/` for their stated
reasons, and that each of the seven audit checks fires on a planted defect.

## Scope of the proofs

The CLR theorems are stated over an abstract commutative ring with a
logarithm-like homomorphism and a reciprocal `ι` of the dimension `n`
(`n · ι = 1`). They are proved about **exact algebra**. The Julia code computes
in `Float64`, so its tests check the same properties **up to a tolerance**. The
proofs say the definitions are right. They do not say the floating-point code
equals the definitions.

## Core layer

| Property | Proved (Agda) | Tested (Julia) | Mutant the test rejects |
|---|---|---|---|
| P1a — CLR components sum to zero | `CLR.clr-sum-zero` | `property_clr.jl` "P1a" | `mutant_clr_uncentred`, `mutant_clr_first_ref` |
| P1b — CLR is scale-invariant | `CLR.clr-scale-invariant` (via `clr-perturb`, `clr-constant`) | "P1b" | `mutant_clr_uncentred` |
| P1c — CLR is permutation-equivariant | `CLR.clr-permute` | "P1c" | `mutant_clr_position_weighted` |
| CLR turns perturbation into addition | `CLR.clr-perturb` | not tested separately | — |
| CLR is **not** invariant under an arbitrary perturbation (n = 2) | `CLR.Instance.clr-not-shift-invariant` | — | — |
| `LogHom` and `Reciprocal` are jointly inhabited for every n ≥ 1, so the CLR theorems are not vacuous (an exact-rational model, not real logarithms) | `CLR.Instance.additive-log`, `CLR.Instance.reciprocal` | — | — |
| BH scaling step is exactly `M·n / ((j+1)(d+1))` | `BenjaminiHochberg.bhScale-numerator`, `bhScale-denominator` | — | — |
| BH q-value is monotone in the family size `M` | `BenjaminiHochberg.monotone-in-family-size` | — | — |
| P5 — BH output ≥ raw, ≤ 1, order-preserving | **not proved** | `property_bh.jl` "P5" | `mutant_bh_no_cummin` |
| BH equals R `p.adjust(, "BH")` | not proved | `unit_core.jl` known answers | — |

P5 and the Agda BH theorems are **different properties**. The Agda file proves
facts about one scaling step; P5 is about the whole adjusted vector, including
the running minimum. Neither substitutes for the other.

## Open obligations

- **BH non-negativity** of the scaled value. Mechanical, not conceptual:
  `+ M * + n` does not reduce to `+ (M ℕ.* n)`, so the stdlib non-negativity
  lemmas do not apply without a transport. It is not asserted anywhere.
- **CLR beyond algebra**: nothing about the real numbers, `exp`, continuity,
  floating point, or the geometry of the simplex is proved.
- **The running minimum** (step-up) of BH is not modelled in Agda.
- **Toolchain**: the proofs need agda-stdlib on the 3.0 development line, pinned
  by SHA in `bootstrap.sh`. They do not compile against v2.4 or v2.1.

## Methods

Not proved, and not claimed. Calibration, consistency and FDR control are
statistical properties, and the plan is to test them against the Bioconductor
reference implementations (ALDEx2, ANCOMBC). That oracle does **not exist yet**:
there is no RCall extension and no `test/oracle_*.jl`. `src/Methods/ALDEx2.jl` is
work in progress and is not yet included in the module.
