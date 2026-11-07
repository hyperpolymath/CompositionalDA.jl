-- SPDX-License-Identifier: MPL-2.0
-- SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
-- EXPECT: 0 != 1
{-# OPTIONS --safe --without-K #-}
module LostAtom where
open import Data.Nat using (_∸_)
open import Relation.Binary.PropositionalEquality using (_≡_; refl)
false-inclusive-lower : 6 ∸ 6 ≡ 1
false-inclusive-lower = refl
