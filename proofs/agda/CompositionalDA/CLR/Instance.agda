-- SPDX-License-Identifier: MPL-2.0
-- SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
--
-- A concrete model of the two CLR assumptions, and a proved counterexample.
--
-- The model is deliberately the simplest one: the ring is the exact rationals
-- `ℚᵘ`, the "positive values" are the naturals under addition, and `log` is the
-- inclusion `a ↦ a/1`.  This is the assumptions *written additively* — it is not
-- a model of real logarithms — and it serves two purposes:
--
--   1. It shows `LogHom` and `Reciprocal` are jointly inhabitable for every
--      `n ≥ 1`, so the theorems of `CompositionalDA.CLR` are not vacuous.
--   2. It hosts `clr-not-shift-invariant`: a *proved* refutation of the
--      tempting false claim that CLR is invariant under perturbation by an
--      arbitrary (non-constant) composition.  In this model perturbation is
--      literally an additive shift of `x`, so this is exactly the "clr is
--      invariant under an additive shift" fallacy, refuted on `n = 2`.
--      (`n = 1` would prove nothing: there every CLR is identically zero.)

{-# OPTIONS --without-K --safe #-}

module CompositionalDA.CLR.Instance where

open import CompositionalDA.Prelude
open import CompositionalDA.LogHom using (LogHom; Reciprocal)
import CompositionalDA.CLR as CLR

open import Algebra.Bundles using (CommutativeRing)
open import Data.Fin.Patterns using (0F; 1F)
open import Data.Integer using (+_)
open import Data.Nat.Base using (ℕ; zero; suc)
open import Data.Nat.Properties using (+-0-commutativeMonoid)
open import Data.Rational.Unnormalised using (ℚᵘ; mkℚᵘ; _≃_)
open import Data.Rational.Unnormalised.Base using (*≡*)
open import Data.Rational.Unnormalised.Properties using (+-*-commutativeRing)
open import Data.Vec.Functional using (Vector)
open import Relation.Nullary.Negation using (¬_)

open import Algebra.Properties.Semiring.Mult
  (CommutativeRing.semiring +-*-commutativeRing) using (_×_)

------------------------------------------------------------------------
-- The additive model of the assumptions

additive-log : LogHom +-*-commutativeRing +-0-commutativeMonoid
additive-log = record
  { log     = λ a → mkℚᵘ (+ a) 0
  ; log-hom = λ a b → ≃-sym (common-denominator-+ (+ a) (+ b) 0)
  }

-- k copies of 1/(m+1) make k/(m+1).
×-unit-fraction : ∀ k m → k × mkℚᵘ (+ 1) m ≃ mkℚᵘ (+ k) m
×-unit-fraction zero    m = ≃-sym (zero-over m)
×-unit-fraction (suc k) m =
  ≃-trans (+-cong-≃ʳ {k × mkℚᵘ (+ 1) m} {mkℚᵘ (+ k) m} {mkℚᵘ (+ 1) m} (×-unit-fraction k m))
          (common-denominator-+ (+ 1) (+ k) m)

-- 1/(m+1) is a reciprocal of m+1.  Every n ≥ 1 is covered; n = 0 has no
-- reciprocal in any ring and no composition with zero parts either.
reciprocal : ∀ m → Reciprocal +-*-commutativeRing (suc m)
reciprocal m = record
  { ι     = mkℚᵘ (+ 1) m
  ; ι-inv = ≃-trans (×-unit-fraction (suc m) m) (whole-is-one m)
  }

------------------------------------------------------------------------
-- The theorems, instantiated at n = 2

module CLR₂ = CLR +-*-commutativeRing additive-log 2 (reciprocal 1)

x₀ : Vector ℕ 2
x₀ _ = 1

y₀ : Vector ℕ 2
y₀ 0F = 1
y₀ 1F = 2

-- Perturbing by a non-constant composition moves the CLR.  This is the
-- statement `proofs/agda/reject/FakeCLR.agda` tries to prove and cannot.
clr-not-shift-invariant : ¬ (CLR₂.clr (CLR₂.perturb x₀ y₀) 0F ≃ CLR₂.clr x₀ 0F)
clr-not-shift-invariant (*≡* ())
