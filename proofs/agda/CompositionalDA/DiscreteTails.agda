-- SPDX-License-Identifier: MPL-2.0
-- SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
{-# OPTIONS --safe --without-K #-}
module CompositionalDA.DiscreteTails where

open import Data.Nat using (ℕ; zero; suc; _+_; _∸_)
open import Relation.Binary.PropositionalEquality using (_≡_; refl; cong)

-- A finite outcome space is partitioned into below, observed, and above.
-- Subtracting the inclusive upper count from the total leaves only below.
cancel-prefix : ∀ upper below → (upper + below) ∸ upper ≡ below
cancel-prefix zero below = refl
cancel-prefix (suc upper) below = cancel-prefix upper below

inclusive-lower : ∀ below observed above →
  (((observed + above) + below) ∸ (observed + above)) + observed ≡
  below + observed
inclusive-lower below observed above =
  cong (_+ observed) (cancel-prefix (observed + above) below)

-- The extreme allocation in a two-versus-two rank test has one observed
-- outcome among six. Complementing its upper tail loses that one outcome.
lost-atom : (6 ∸ 6) + 1 ≡ 1
lost-atom = refl

strict-lower : 6 ∸ 6 ≡ 0
strict-lower = refl
