-- SPDX-License-Identifier: MPL-2.0
-- SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
{-# OPTIONS --safe --without-K #-}

-- The finite aggregation inequality behind ALDEx2's consistent-sign rule.
-- This is algebra over exact rationals, not a calibration/FDR theorem.
module CompositionalDA.PosteriorP where

open import Data.Nat.Base using (ℕ; zero; suc)
open import Data.Vec.Base using (Vec; []; _∷_; zipWith)
open import Data.Rational.Unnormalised.Base using (ℚᵘ; 0ℚᵘ; _+_; _⊓_; _≤_)
open import Data.Rational.Unnormalised.Properties
  using (≤-refl; +-mono-≤; p⊓q≤p; p⊓q≤q; ⊓-glb)

sum : ∀ {n} → Vec ℚᵘ n → ℚᵘ
sum [] = 0ℚᵘ
sum (x ∷ xs) = x + sum xs

sum-min≤left : ∀ {n} (xs ys : Vec ℚᵘ n) →
  sum (zipWith _⊓_ xs ys) ≤ sum xs
sum-min≤left [] [] = ≤-refl
sum-min≤left (x ∷ xs) (y ∷ ys) =
  +-mono-≤ (p⊓q≤p x y) (sum-min≤left xs ys)

sum-min≤right : ∀ {n} (xs ys : Vec ℚᵘ n) →
  sum (zipWith _⊓_ xs ys) ≤ sum ys
sum-min≤right [] [] = ≤-refl
sum-min≤right (x ∷ xs) (y ∷ ys) =
  +-mono-≤ (p⊓q≤q x y) (sum-min≤right xs ys)

-- Multiplication by the common positive 1/mc gives the corresponding
-- inequality of means. This theorem does not assume either vector is random.
consistent-sign-bound : ∀ {n} (xs ys : Vec ℚᵘ n) →
  sum (zipWith _⊓_ xs ys) ≤ (sum xs ⊓ sum ys)
consistent-sign-bound xs ys =
  ⊓-glb (sum-min≤left xs ys) (sum-min≤right xs ys)
