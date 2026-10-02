-- SPDX-License-Identifier: MPL-2.0
-- SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
--
-- The two named assumptions under which the CLR theorems are proved.
--
-- Neither is a postulate.  Each is a record whose fields a caller must supply,
-- so every theorem in `CompositionalDA.CLR` is a function of a witness and
-- nothing is asserted without proof.  `CompositionalDA.CLR.Instance` shows the
-- two records are jointly inhabitable, so the theorems are not vacuous.
--
-- `LogHom` is ported from the MetaManifold proof tree
-- (`MetaManifold.ILR.Invariance`): a logarithm is any map from a commutative
-- monoid of "positive" values into the ring that turns the monoid operation
-- into ring addition.  That is the only property of `log` any theorem here
-- uses; nothing about real numbers, continuity or `exp` is assumed.
--
-- `Reciprocal` is new to this tree: CLR divides by the number of parts `n`, so
-- the ring must contain an element `ι` with `n · ι = 1`.  Over the rationals
-- `ι = 1/n` for `n ≥ 1`; over the integers no such element exists for `n ≥ 2`,
-- which is exactly why `proofs/agda/reject/MissingReciprocal.agda` fails.

{-# OPTIONS --without-K --safe #-}

module CompositionalDA.LogHom where

open import Algebra.Bundles using (CommutativeRing; CommutativeMonoid)
open import Data.Nat.Base using (ℕ)
open import Level using (_⊔_)

------------------------------------------------------------------------
-- The log seam

record LogHom {c ℓ m ℓm} (R : CommutativeRing c ℓ) (M : CommutativeMonoid m ℓm)
       : Set (c ⊔ ℓ ⊔ m ⊔ ℓm) where
  open CommutativeRing R using (Carrier; _≈_; _+_)
  open CommutativeMonoid M using () renaming (Carrier to Pos; _∙_ to _×ₘ_)
  field
    log     : Pos → Carrier
    log-hom : ∀ a b → log (a ×ₘ b) ≈ log a + log b

------------------------------------------------------------------------
-- A reciprocal of the number of parts

module _ {c ℓ} (R : CommutativeRing c ℓ) where
  open CommutativeRing R using (Carrier; _≈_; 1#; semiring)
  open import Algebra.Properties.Semiring.Mult semiring using (_×_)

  record Reciprocal (n : ℕ) : Set (c ⊔ ℓ) where
    field
      ι     : Carrier
      ι-inv : n × ι ≈ 1#
