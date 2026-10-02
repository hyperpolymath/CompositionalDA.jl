-- SPDX-License-Identifier: MPL-2.0
-- SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
--
-- Negative control: a CLR "theorem" that is false must not type-check.
--
-- The claim: CLR is invariant under perturbation of x by an arbitrary
-- composition y — the "clr is invariant under an additive shift of x"
-- fallacy.  The true theorem (`clr-scale-invariant`) needs y to be constant;
-- this drops that hypothesis.  `CompositionalDA.CLR.Instance` proves the
-- negation of exactly this statement at n = 2, so no proof can exist; the
-- attempt below is the one a careless author would try, and Agda reports the
-- two cross-multiplied integers that would have to be equal (`-2` and `0`).
--
-- EXPECT: negsuc [0-9]+\) != \(Agda\.Builtin\.Int\.Int\.pos 0\)

{-# OPTIONS --without-K --safe #-}

module reject.FakeCLR where

open import CompositionalDA.CLR.Instance using (module CLR₂; x₀; y₀)
open import Data.Fin.Patterns using (0F)
open import Data.Rational.Unnormalised using (_≃_)
open import Data.Rational.Unnormalised.Base using (*≡*)
open import Relation.Binary.PropositionalEquality using (refl)

clr-shift-invariant : CLR₂.clr (CLR₂.perturb x₀ y₀) 0F ≃ CLR₂.clr x₀ 0F
clr-shift-invariant = *≡* refl
