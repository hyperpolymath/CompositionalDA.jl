-- SPDX-License-Identifier: MPL-2.0
-- SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
--
-- P1: the centred log-ratio (CLR) transform, over the named assumptions in
-- `CompositionalDA.LogHom`.
--
-- A composition with `n` parts is a vector `Fin n → Pos` of "positive" values
-- in a commutative monoid `M` (the monoid operation `×ₘ` is Aitchison's
-- perturbation, componentwise).  A logarithm is any `LogHom`: a map into a
-- commutative ring `R` that sends `×ₘ` to `+`.  Dividing by `n` needs a
-- reciprocal `ι` with `n · ι = 1`; that is the second record.
--
--   clr x i = log (x i) − ι · Σⱼ log (x j)
--
-- Proved here, each as a function of the two witnesses and nothing else:
--
--   (a) clr-sum-zero          Σᵢ clr x i = 0
--   (b) clr-scale-invariant   clr (c · x) = clr x          for every c : Pos
--   (c) clr-permute           clr (x ∘ σ) = clr x ∘ σ      for every permutation σ
--
-- plus the linearity `clr-perturb` (clr turns perturbation into addition) and
-- `clr-constant` (a constant composition has zero CLR) that (b) is built from.
--
-- What is NOT here: anything about real numbers, `exp`, continuity, floating
-- point, or the geometry of the simplex.  See `proofs/PROOF-STATUS.md`.

{-# OPTIONS --without-K --safe #-}

open import Algebra.Bundles using (CommutativeRing; CommutativeMonoid)
open import Data.Nat.Base using (ℕ)
open import CompositionalDA.LogHom using (LogHom; Reciprocal)

module CompositionalDA.CLR
  {c ℓ m ℓm} (R : CommutativeRing c ℓ) {M : CommutativeMonoid m ℓm}
  (L : LogHom R M) (n : ℕ) (ρ : Reciprocal R n)
  where

open import Data.Fin.Base using (Fin)
open import Data.Fin.Permutation using (Permutation; _⟨$⟩ʳ_)
open import Data.Maybe.Base using (nothing)
open import Data.Vec.Functional using (Vector; map; replicate; rearrange; zipWith)

open CommutativeRing R
open CommutativeMonoid M using () renaming (Carrier to Pos; _∙_ to _×ₘ_)
open LogHom L
open Reciprocal ρ

open import Algebra.Properties.CommutativeMonoid.Sum +-commutativeMonoid
  using (sum; ∑-distrib-+; sum-replicate; sum-permute; sum-cong-≋)
open import Algebra.Properties.Semiring.Mult semiring
  using (_×_; ×-congʳ; ×-assoc-*; ×-comm-*)
open import Algebra.Properties.Ring ring using (-‿distribʳ-*)
open import Relation.Binary.Reasoning.Setoid setoid

import Tactic.RingSolver.Core.AlmostCommutativeRing as ACR
open import Tactic.RingSolver.NonReflective (ACR.fromCommutativeRing R (λ _ → nothing))
  using (solve; _⊜_; _⊕_; _⊗_; ⊝_)

------------------------------------------------------------------------
-- Definitions

-- Centring: subtract ι times the total from every entry.
centre : Vector Carrier n → Vector Carrier n
centre v i = v i - ι * sum v

-- The centred log-ratio transform.
clr : Vector Pos n → Vector Carrier n
clr x = centre (map log x)

-- Perturbation, the compositional "addition": componentwise ×ₘ.
perturb : Vector Pos n → Vector Pos n → Vector Pos n
perturb = zipWith _×ₘ_

-- Scaling by a constant c : Pos is perturbation by the constant composition.
scale : Pos → Vector Pos n → Vector Pos n
scale k x = perturb (replicate n k) x

------------------------------------------------------------------------
-- (a) The CLR components sum to zero

centre-sum-zero : ∀ v → sum (centre v) ≈ 0#
centre-sum-zero v = begin
  sum (centre v)                         ≈⟨ ∑-distrib-+ v (λ _ → - (ι * s)) ⟩
  sum v + sum (replicate n (- (ι * s)))  ≈⟨ +-congˡ (sum-replicate n) ⟩
  sum v + n × (- (ι * s))                ≈⟨ +-congˡ (×-congʳ n (-‿distribʳ-* ι s)) ⟩
  sum v + n × (ι * - s)                  ≈⟨ +-congˡ (×-assoc-* n ι (- s)) ⟨
  sum v + (n × ι) * - s                  ≈⟨ +-congˡ (*-congʳ ι-inv) ⟩
  sum v + 1# * - s                       ≈⟨ +-congˡ (*-identityˡ (- s)) ⟩
  sum v + - sum v                        ≈⟨ -‿inverseʳ (sum v) ⟩
  0#                                     ∎
  where s = sum v

clr-sum-zero : ∀ (x : Vector Pos n) → sum (clr x) ≈ 0#
clr-sum-zero x = centre-sum-zero (map log x)

------------------------------------------------------------------------
-- Linearity: clr turns perturbation into addition

clr-perturb : ∀ (x y : Vector Pos n) i →
              clr (perturb x y) i ≈ clr x i + clr y i
clr-perturb x y i = begin
  log (x i ×ₘ y i) - ι * sum (map log (perturb x y))
    ≈⟨ +-cong (log-hom (x i) (y i)) (-‿cong (*-congˡ sum-logs)) ⟩
  (a + b) - ι * (s + t)
    ≈⟨ solve 5 (λ a b s t ι → ((a ⊕ b) ⊕ ⊝ (ι ⊗ (s ⊕ t)))
                              ⊜ ((a ⊕ ⊝ (ι ⊗ s)) ⊕ (b ⊕ ⊝ (ι ⊗ t))))
             refl a b s t ι ⟩
  (a - ι * s) + (b - ι * t)
    ∎
  where
  a = log (x i)
  b = log (y i)
  s = sum (map log x)
  t = sum (map log y)
  sum-logs : sum (map log (perturb x y)) ≈ s + t
  sum-logs = trans (sum-cong-≋ (λ j → log-hom (x j) (y j)))
                   (∑-distrib-+ (map log x) (map log y))

-- A constant composition has zero CLR.
clr-constant : ∀ (k : Pos) i → clr (replicate n k) i ≈ 0#
clr-constant k i = begin
  log k - ι * sum (replicate n (log k))  ≈⟨ +-congˡ (-‿cong (*-congˡ (sum-replicate n))) ⟩
  log k - ι * (n × log k)                ≈⟨ +-congˡ (-‿cong (×-comm-* n ι (log k))) ⟩
  log k - n × (ι * log k)                ≈⟨ +-congˡ (-‿cong (×-assoc-* n ι (log k))) ⟨
  log k - (n × ι) * log k                ≈⟨ +-congˡ (-‿cong (*-congʳ ι-inv)) ⟩
  log k - 1# * log k                     ≈⟨ +-congˡ (-‿cong (*-identityˡ (log k))) ⟩
  log k - log k                          ≈⟨ -‿inverseʳ (log k) ⟩
  0#                                     ∎

------------------------------------------------------------------------
-- (b) Scale invariance: clr (c · x) = clr x

clr-scale-invariant : ∀ (k : Pos) (x : Vector Pos n) i → clr (scale k x) i ≈ clr x i
clr-scale-invariant k x i = begin
  clr (perturb (replicate n k) x) i   ≈⟨ clr-perturb (replicate n k) x i ⟩
  clr (replicate n k) i + clr x i     ≈⟨ +-congʳ (clr-constant k i) ⟩
  0# + clr x i                        ≈⟨ +-identityˡ (clr x i) ⟩
  clr x i                             ∎

------------------------------------------------------------------------
-- (c) Permutation equivariance: clr (x ∘ σ) = clr x ∘ σ

clr-permute : ∀ (π : Permutation n n) (x : Vector Pos n) i →
              clr (rearrange (π ⟨$⟩ʳ_) x) i ≈ clr x (π ⟨$⟩ʳ i)
clr-permute π x i = +-congˡ (-‿cong (*-congˡ (sym (sum-permute (map log x) π))))
