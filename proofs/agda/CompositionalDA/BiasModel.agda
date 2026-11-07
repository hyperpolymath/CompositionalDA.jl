-- SPDX-License-Identifier: MPL-2.0
-- SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
{-# OPTIONS --safe --without-K #-}

open import Algebra.Bundles using (CommutativeRing)

module CompositionalDA.BiasModel {c ℓ} (R : CommutativeRing c ℓ) where

open CommutativeRing R
open import Relation.Binary.Reasoning.Setoid setoid

-- The observable linear predictor cannot distinguish an abundance shift
-- from an opposite sampling-fraction shift. An identification assumption is
-- necessary; count data alone cannot discharge it.
observational-shift : ∀ μ d s → (μ + s) + (d - s) ≈ μ + d
observational-shift μ d s = begin
  (μ + s) + (d - s)
    ≈⟨ +-cong refl (+-comm d (- s)) ⟩
  (μ + s) + ((- s) + d)
    ≈⟨ sym (+-assoc (μ + s) (- s) d) ⟩
  ((μ + s) + (- s)) + d
    ≈⟨ +-cong (+-assoc μ s (- s)) refl ⟩
  (μ + (s + (- s))) + d
    ≈⟨ +-cong (+-cong refl (-‿inverseʳ s)) refl ⟩
  (μ + 0#) + d
    ≈⟨ +-cong (+-identityʳ μ) refl ⟩
  μ + d ∎

-- Subtraction removes EXACTLY the supplied bias. No claim that the estimated
-- bias equals the population bias is assumed to have been proved here.
bias-cancellation : ∀ μ ν δ → ((μ - ν) + δ) - δ ≈ μ - ν
bias-cancellation μ ν δ = begin
  ((μ - ν) + δ) - δ
    ≈⟨ +-assoc (μ - ν) δ (- δ) ⟩
  (μ - ν) + (δ + (- δ))
    ≈⟨ +-cong refl (-‿inverseʳ δ) ⟩
  (μ - ν) + 0#
    ≈⟨ +-identityʳ (μ - ν) ⟩
  μ - ν ∎

-- Residual bias is precisely estimation error, with the opposite sign.
residual-bias : ∀ μ ν δ estimate →
  ((μ - ν) + δ) - estimate ≈ (μ - ν) + (δ - estimate)
residual-bias μ ν δ estimate = +-assoc (μ - ν) δ (- estimate)
