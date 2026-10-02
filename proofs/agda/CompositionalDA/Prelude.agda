-- SPDX-License-Identifier: MPL-2.0
-- SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
--
-- Exact-rational foundations for CompositionalDA's proof layer.
--
-- Ported from the MetaManifold proof tree (`MetaManifold.Prelude`), reduced to
-- the rational core that P1 (CLR) and P5 (Benjamini–Hochberg) use.  The
-- refusal/outcome types and list sums of the original are not needed here and
-- are deliberately not carried over.
--
-- `ℚᵘ` is the *unnormalised* rationals: an integer numerator and a positive
-- denominator, with equality `_≃_` by cross-multiplication.  Every theorem is
-- stated about values (`_≃_`), so each transfers unchanged to the normalised
-- `Rational{BigInt}` the Julia layer holds.  The denominator is `suc _` by
-- construction, so a zero denominator is unrepresentable rather than checked.
--
-- There are no `postulate`s anywhere in `proofs/agda/CompositionalDA`;
-- `proofs/tests/axiom-audit.sh` checks that, and checks its own control.

{-# OPTIONS --without-K --safe #-}

module CompositionalDA.Prelude where

open import Data.Integer as ℤ using (ℤ; +_; +0; +[1+_]; -[1+_]; _*_; _+_; _≤_)
open import Data.Integer.Properties as ℤₚ
  using (pos-*; pos-+; *-comm; *-zeroˡ; _≡?_; *-monoʳ-≤-nonNeg)
  renaming (+-*-commutativeRing to +-*-commutativeRingℤ)
open import Data.Nat.Base using (z≤n)
open import Data.Nat.Base as ℕ using (ℕ; zero; suc)
open import Data.Rational.Unnormalised
  using (ℚᵘ; mkℚᵘ; _≃_; ↥_; ↧_; ↧ₙ_; 0ℚᵘ; 1ℚᵘ)
  renaming (_+_ to _+ℚ_; _*_ to _*ℚ_; _/_ to _/ℚ_; _≤_ to _≤ℚ_)
open import Data.Rational.Unnormalised.Base using (*≡*; *≤*)
open import Data.Rational.Unnormalised.Properties
  using (+-*-commutativeRing; ≃-isEquivalence; +-isMagma; *-isMagma)
open import Relation.Binary using (Decidable)
open import Relation.Binary.PropositionalEquality
  using (_≡_; _≢_; refl; sym; trans; cong)
open import Relation.Nullary.Decidable using (yes; no)

------------------------------------------------------------------------
-- Exact rationals

infix 4 _≟ℚᵘ_

_≟ℚᵘ_ : Decidable _≃_
p ≟ℚᵘ q with (↥ p * ↧ q) ≡? (↥ q * ↧ p)
... | yes e  = yes (*≡* e)
... | no  ¬e = no (λ where (*≡* e) → ¬e e)

-- Setoid structure over `_≃_`.
open import Algebra.Structures using (IsMagma)
open import Relation.Binary.Structures using (IsEquivalence)
open IsEquivalence ≃-isEquivalence public
  using () renaming (refl to ≃-refl; sym to ≃-sym; trans to ≃-trans)

-- Congruence of the two operations, taken from the library's own structures
-- rather than re-proved.
open IsMagma +-isMagma public using () renaming (∙-cong to +-cong-≃)
open IsMagma *-isMagma public using () renaming (∙-cong to *-cong-≃)

-- Ring solvers, so that rational and integer identities are discharged by
-- computation rather than by hand-written chains.
open import Algebra.Solver.Ring.AlmostCommutativeRing using (fromCommutativeRing)
import Algebra.Solver.Ring.Simple as RingSolver

module ℚᵘ-Solver = RingSolver (fromCommutativeRing +-*-commutativeRing) _≟ℚᵘ_
open ℚᵘ-Solver public using (con; _:+_; _:*_; _:=_) renaming (solve to solve-ℚ)

module ℤ-Solver = RingSolver (fromCommutativeRing +-*-commutativeRingℤ) _≡?_
open ℤ-Solver public using ()
  renaming (solve to solve-ℤ; _:=_ to _:=ℤ_; con to conℤ; _:+_ to _:+ℤ_; _:*_ to _:*ℤ_)

-- `n/d` with the positivity proof made explicit.
_/[_]_ : (n : ℤ) (d : ℕ) → d ≢ 0 → ℚᵘ
n /[ d ] p = _/ℚ_ n d {{ℕ.≢-nonZero p}}

infixl 7 _/[_]_

-- `a/d + b/d ≃ (a + b)/d`: quotients over a common denominator add by adding
-- their numerators.
common-denominator-+ : ∀ (a b : ℤ) (d : ℕ) →
                       mkℚᵘ a d +ℚ mkℚᵘ b d ≃ mkℚᵘ (a + b) d
common-denominator-+ a b d = *≡* cross
  where
  -- ↥(p+q) · ↧r  ≡  ↥r · ↧(p+q), where p = a/d, q = b/d, r = (a+b)/d:
  --   (a·(1+d) + b·(1+d)) · (1+d)  ≡  (a+b) · ((1+d)·(1+d))
  cross : (a * +[1+ d ] + b * +[1+ d ]) * +[1+ d ]
          ≡ (a + b) * +[1+ d ℕ.+ d ℕ.* suc d ]
  cross = trans
    (solve-ℤ 3 (λ a b d → (a :*ℤ d :+ℤ b :*ℤ d) :*ℤ d :=ℤ (a :+ℤ b) :*ℤ (d :*ℤ d))
       refl a b (+[1+ d ]))
    (cong ((a + b) *_) (sym (ℤₚ.pos-* (suc d) (suc d))))

-- `n/n ≃ 1`.  No side condition that `n` is non-zero: in `ℚᵘ` the denominator
-- is `suc _`, so `n/n` here means `(n+1)/(n+1)` and the zero-denominator case
-- is not a case this type has.
whole-is-one : ∀ (n : ℕ) → mkℚᵘ (+[1+ n ]) n ≃ 1ℚᵘ
whole-is-one n = *≡* (*-comm (+[1+ n ]) (+[1+ 0 ]))

-- `0/n ≃ 0`.
zero-over : ∀ (d : ℕ) → mkℚᵘ (+ 0) d ≃ 0ℚᵘ
zero-over d = *≡* (trans (ℤₚ.*-zeroˡ (+[1+ 0 ])) (sym (ℤₚ.*-zeroˡ (+[1+ d ]))))

-- Right-congruence of `+ℚ`, with the implicit arguments in the order callers
-- find convenient.
+-cong-≃ʳ : ∀ {u v w : ℚᵘ} → u ≃ v → (w +ℚ u) ≃ (w +ℚ v)
+-cong-≃ʳ {u} {v} {w} p = +-cong-≃ {w} {w} {u} {v} ≃-refl p

-- Propositional equality implies value equality.
≡⇒≃ : ∀ {p q : ℚᵘ} → p ≡ q → p ≃ q
≡⇒≃ refl = ≃-refl

-- Over a fixed positive denominator, the order on rationals is the order on
-- numerators.  This is the whole content of P5's monotonicity result.
numerator-monotone-≤ : ∀ (n m : ℤ) (d : ℕ) → n ≤ m → mkℚᵘ n d ≤ℚ mkℚᵘ m d
numerator-monotone-≤ n m d n≤m =
  *≤* (ℤₚ.*-monoʳ-≤-nonNeg +[1+ d ] {{ℤ.nonNegative (ℤ.+≤+ z≤n)}} n≤m)
