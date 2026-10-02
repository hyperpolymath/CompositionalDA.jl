-- SPDX-License-Identifier: MPL-2.0
-- SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
--
-- Negative control: the `Reciprocal` hypothesis cannot be waved away.
--
-- `CompositionalDA.CLR` needs an element ι with `n · ι = 1`.  Over the
-- integers with n = 2 there is no such ι, so the only way to instantiate the
-- theorems is to fake the witness.  Here the fake is `ι = 1` with the law
-- discharged by `refl`, which would need `2 · 1 = 1` to hold definitionally in
-- ℤ.  It does not: `2 · 1` normalises to `+ 2`, and Agda reports the mismatch
-- on the underlying naturals, `2 != 1`.
--
-- EXPECT: ^2 != 1 of type Agda\.Builtin\.Nat\.Nat

{-# OPTIONS --without-K --safe #-}

module reject.MissingReciprocal where

open import CompositionalDA.LogHom using (Reciprocal)
open import Data.Integer using (+_)
open import Data.Integer.Properties using (+-*-commutativeRing)
open import Relation.Binary.PropositionalEquality using (refl)

no-half : Reciprocal +-*-commutativeRing 2
no-half = record { ι = + 1 ; ι-inv = refl }
