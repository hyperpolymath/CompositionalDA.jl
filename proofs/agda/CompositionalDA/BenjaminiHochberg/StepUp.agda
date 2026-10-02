-- SPDX-License-Identifier: MPL-2.0
-- SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
--
-- Benjamini–Hochberg step-up, as the running minimum that `src/Core/bh.jl`
-- computes (P5).
--
-- The Julia loop, after sorting the p-values into DESCENDING order, is
--
--     running = 1.0
--     for (r, k) in enumerate(ord)          # r = 1 is the largest p-value
--         rank    = m - r + 1                # m, m-1, …, 1
--         running = min(running, q[k] * m / rank)
--         adj[k]  = running
--     end
--
-- `stepUp` below is that loop, on exact rationals.  The rank of the head is
-- the number of p-values still to process, so it counts down `m, m-1, …, 1`
-- exactly as `m - r + 1` does, and `running` starts at `1ℚᵘ`.  `adjust` fixes
-- the family size to the length of the vector, which is what `m = length(idx)`
-- does after the NaNs are dropped.
--
-- What is proved about `adjust`, each as a statement about `lookup` at a
-- position of the sorted sequence, and then re-indexed through the sorting
-- permutation so that it reads at the ORIGINAL index the way the Julia
-- property test `test/property_bh.jl` reads it:
--
--   * `adjust-≤1`        — every adjusted value is ≤ 1.          (no hypothesis)
--   * `adjust-descending`— the output is descending too.         (no hypothesis)
--   * `adjust-≥`         — q_i ≥ p_i.      Needs: input descending, 0 ≤ p ≤ 1.
--   * `adjust-monotone`  — p_i ≤ p_j ⇒ q_i ≤ q_j, at ANY two positions,
--                          ties included.  Needs: input descending, 0 ≤ p.
--
-- The hypotheses are not formalities.  `0 ≤ p ≤ 1` is the `ArgumentError`
-- guard in `bh_adjust`; without `0 ≤ p` the scaling `p · m/rank` runs the
-- wrong way, and without `p ≤ 1` the initial `running = 1` already undercuts
-- the largest p-value.  "Input descending" is what `sortperm(q; rev = true)`
-- establishes; it is assumed here, not re-proved, and that assumption is
-- recorded in `proofs/PROOF-STATUS.md`.
--
-- Ties are the part a test on distinct values would never exercise: equal
-- p-values occupy consecutive ranks and so receive DIFFERENT scaled values,
-- and the proof that they nevertheless receive the same adjusted value
-- (`stepUp-tie`) is where `0 ≤ p` is actually used.
--
-- Everything is over `ℚᵘ`, the unnormalised rationals, so no gcd step and no
-- rounding sits between the arithmetic and the proof.  `_⊓_` on `ℚᵘ` does
-- not reduce on variables, so every step below goes through the lattice
-- lemmas of `Data.Rational.Unnormalised.Properties`, never through `refl`.

{-# OPTIONS --without-K --safe #-}

module CompositionalDA.BenjaminiHochberg.StepUp where

open import Data.Nat.Base as ℕ using (ℕ; zero; suc; z≤n; s≤s)
open import Data.Nat.Properties as ℕₚ using (n≤1+n; <⇒≤)
open import Data.Integer.Base as ℤ using (ℤ; +_; +≤+)
open import Data.Integer.Properties as ℤₚ using ()
open import Data.Fin.Base using (Fin; zero; suc; toℕ)
open import Data.Fin.Permutation using (Permutation; _⟨$⟩ʳ_; _⟨$⟩ˡ_; inverseʳ)
open import Data.Vec.Base using (Vec; []; _∷_; lookup; tabulate)
open import Data.Vec.Properties using (lookup∘tabulate)
open import Data.Vec.Relation.Unary.All as All using (All; []; _∷_)
open import Data.Vec.Relation.Unary.All.Properties using (tabulate⁺)
open import Data.Vec.Relation.Unary.Linked using (Linked; []; [-]; _∷_)
open import Data.Sum.Base using (inj₁; inj₂)
open import Relation.Binary.PropositionalEquality
  using (_≡_; refl; sym; trans; cong; subst; subst₂)
open import Data.Rational.Unnormalised.Base
  using (ℚᵘ; mkℚᵘ; 0ℚᵘ; 1ℚᵘ; _≃_; _≤_; _≥_; _*_; _⊓_; *≤*; nonNegative)
open import Data.Rational.Unnormalised.Properties
  using ( ≃-refl; ≃-sym; ≃-trans
        ; ≤-refl; ≤-reflexive; ≤-trans; ≤-antisym
        ; *-identityʳ; *-congʳ; *-monoʳ-≤-nonNeg
        ; p⊓q≤p; p⊓q≤q; ⊓-glb; p≤q⇒p⊓q≃p )

------------------------------------------------------------------------
-- The model

-- `p · m/(k+1)`: the Julia `q[k] * m / rank` with `rank = k + 1`.  `ℚᵘ`
-- stores a denominator-minus-one, so the fraction `m/(k+1)` is `mkℚᵘ (+ m) k`.
scale : ℕ → ℕ → ℚᵘ → ℚᵘ
scale m k p = p * mkℚᵘ (+ m) k

-- The loop.  `run` is Julia's `running`.  The head of a vector of length
-- `suc k` has rank `suc k`, i.e. `scale m k`.
stepUp : (m : ℕ) → ℚᵘ → ∀ {k} → Vec ℚᵘ k → Vec ℚᵘ k
stepUp m run []               = []
stepUp m run {suc k} (p ∷ ps) =
  (run ⊓ scale m k p) ∷ stepUp m (run ⊓ scale m k p) ps

-- `bh_adjust` on the sorted, NaN-free p-values: family size = length,
-- `running` starts at 1.
adjust : ∀ {m} → Vec ℚᵘ m → Vec ℚᵘ m
adjust {m} = stepUp m 1ℚᵘ

------------------------------------------------------------------------
-- The scaling factor

-- `m/(k+1) ≥ 1` whenever the rank `k+1` does not exceed the family size.
1≤ratio : ∀ m k → suc k ℕ.≤ m → 1ℚᵘ ≤ mkℚᵘ (+ m) k
1≤ratio m k sk≤m = *≤*
  (ℤₚ.≤-trans (ℤₚ.≤-reflexive (ℤₚ.*-identityˡ (+ suc k)))
  (ℤₚ.≤-trans (+≤+ sk≤m)
              (ℤₚ.≤-reflexive (sym (ℤₚ.*-identityʳ (+ m))))))

-- A smaller rank gives a larger factor: `m/(k+2) ≤ m/(k+1)`.
ratio-antitone : ∀ m k → mkℚᵘ (+ m) (suc k) ≤ mkℚᵘ (+ m) k
ratio-antitone m k = *≤* (ℤₚ.*-monoˡ-≤-nonNeg (+ m) (+≤+ (n≤1+n (suc k))))

-- Scaling a non-negative p-value by a factor ≥ 1 cannot shrink it.
scale-≥ : ∀ m k p → suc k ℕ.≤ m → 0ℚᵘ ≤ p → p ≤ scale m k p
scale-≥ m k p sk≤m 0≤p =
  ≤-trans (≤-reflexive (≃-sym (*-identityʳ p)))
          (*-monoʳ-≤-nonNeg p {{nonNegative 0≤p}} (1≤ratio m k sk≤m))

-- Moving a non-negative p-value one rank down (towards rank 1) scales it up.
scale-antitone : ∀ m k p → 0ℚᵘ ≤ p → scale m (suc k) p ≤ scale m k p
scale-antitone m k p 0≤p =
  *-monoʳ-≤-nonNeg p {{nonNegative 0≤p}} (ratio-antitone m k)

------------------------------------------------------------------------
-- The running minimum, position by position

-- The output never rises: each value is the previous one `⊓` something.
stepUp-descending : ∀ m run {k} (ps : Vec ℚᵘ k) → Linked _≥_ (stepUp m run ps)
stepUp-descending m run []            = []
stepUp-descending m run (p ∷ [])      = [-]
stepUp-descending m run {suc (suc k)} (p ∷ p′ ∷ ps) =
  p⊓q≤p (run ⊓ scale m (suc k) p) (scale m k p′)
    ∷ stepUp-descending m (run ⊓ scale m (suc k) p) (p′ ∷ ps)

-- Nothing exceeds the initial `running`.
stepUp-≤run : ∀ m run {k} (ps : Vec ℚᵘ k) (i : Fin k) →
              lookup (stepUp m run ps) i ≤ run
stepUp-≤run m run {suc k} (p ∷ ps) zero    = p⊓q≤p run (scale m k p)
stepUp-≤run m run {suc k} (p ∷ ps) (suc i) =
  ≤-trans (stepUp-≤run m (run ⊓ scale m k p) ps i) (p⊓q≤p run (scale m k p))

-- In a descending vector a later position holds a smaller-or-equal value.
desc-lookup : ∀ {k} {xs : Vec ℚᵘ k} → Linked _≥_ xs →
              ∀ (i j : Fin k) → toℕ i ℕ.≤ toℕ j → lookup xs j ≤ lookup xs i
desc-lookup lnk zero zero _ = ≤-refl
desc-lookup {xs = x ∷ x′ ∷ xs} (x′≤x ∷ lnk) zero (suc j) _ =
  ≤-trans (desc-lookup lnk zero j z≤n) x′≤x
desc-lookup {xs = x ∷ xs} (_ ∷ lnk) (suc i) (suc j) (s≤s i≤j) =
  desc-lookup lnk i j i≤j

-- q_i ≥ p_i, given: the rank of the head does not exceed the family size,
-- the head is below the incoming `running`, the input is descending and
-- non-negative.  The last three are what the recursion has to re-establish
-- for the tail: the new `running` is ≥ the head, which is ≥ the next value.
stepUp-≥ : ∀ m run {k} (p : ℚᵘ) (ps : Vec ℚᵘ k) → suc k ℕ.≤ m → p ≤ run →
           Linked _≥_ (p ∷ ps) → All (0ℚᵘ ≤_) (p ∷ ps) →
           ∀ (i : Fin (suc k)) → lookup (p ∷ ps) i ≤ lookup (stepUp m run (p ∷ ps)) i
stepUp-≥ m run {k} p ps sk≤m p≤run _ (0≤p ∷ _) zero =
  ⊓-glb p≤run (scale-≥ m k p sk≤m 0≤p)
stepUp-≥ m run {suc k} p (p′ ∷ ps) sk≤m p≤run (p′≤p ∷ lnk) (0≤p ∷ all) (suc i) =
  stepUp-≥ m (run ⊓ scale m (suc k) p) p′ ps (<⇒≤ sk≤m)
           (≤-trans p′≤p (⊓-glb p≤run (scale-≥ m (suc k) p sk≤m 0≤p)))
           lnk all i

-- Ties.  If the value at position `j` is ≥ the head, then (the input being
-- descending) every value from the head to `j` equals the head, and the
-- running minimum does not move across that block: each later scaled value
-- is ≥ the head's scaled value (`scale-antitone`), which is already ≥ the
-- minimum.  So the adjusted value at `j` IS the adjusted value at the head.
stepUp-tie : ∀ m run {k} (p : ℚᵘ) (ps : Vec ℚᵘ k) →
             Linked _≥_ (p ∷ ps) → All (0ℚᵘ ≤_) (p ∷ ps) →
             ∀ (j : Fin (suc k)) → p ≤ lookup (p ∷ ps) j →
             lookup (stepUp m run (p ∷ ps)) j ≃ run ⊓ scale m k p
stepUp-tie m run p ps _ _ zero _ = ≃-refl
stepUp-tie m run {suc k} p (p′ ∷ ps) (p′≤p ∷ lnk) (0≤p ∷ 0≤p′ ∷ all) (suc j) p≤pⱼ =
  ≃-trans (stepUp-tie m q p′ ps lnk (0≤p′ ∷ all) j (≤-trans p′≤p p≤pⱼ))
          (p≤q⇒p⊓q≃p q≤s′)
  where
    q : ℚᵘ
    q = run ⊓ scale m (suc k) p
    pⱼ≤p′ : lookup (p′ ∷ ps) j ≤ p′
    pⱼ≤p′ = desc-lookup lnk zero j z≤n
    p′≃p : p′ ≃ p
    p′≃p = ≤-antisym p′≤p (≤-trans p≤pⱼ pⱼ≤p′)
    q≤s′ : q ≤ scale m k p′
    q≤s′ = ≤-trans (p⊓q≤q run (scale m (suc k) p))
           (≤-trans (scale-antitone m k p 0≤p)
                    (≤-reflexive (*-congʳ (≃-sym p′≃p))))

-- The tie lemma, started at any position `i` rather than at the head.
stepUp-tie-from : ∀ m run {k} (ps : Vec ℚᵘ k) → Linked _≥_ ps → All (0ℚᵘ ≤_) ps →
                  ∀ (i j : Fin k) → toℕ i ℕ.≤ toℕ j → lookup ps i ≤ lookup ps j →
                  lookup (stepUp m run ps) i ≃ lookup (stepUp m run ps) j
stepUp-tie-from m run (p ∷ ps) lnk all zero j _ p≤pⱼ =
  ≃-sym (stepUp-tie m run p ps lnk all j p≤pⱼ)
stepUp-tie-from m run {suc k} (p ∷ ps) (_ ∷ lnk) (_ ∷ all) (suc i) (suc j) (s≤s i≤j) pᵢ≤pⱼ =
  stepUp-tie-from m (run ⊓ scale m k p) ps lnk all i j i≤j pᵢ≤pⱼ

-- p_i ≤ p_j ⇒ q_i ≤ q_j at any two positions.  Either `j` is at or before
-- `i`, and the output is descending; or `i` is at or before `j`, in which
-- case the input (also descending) is constant between them and the tie
-- lemma makes the outputs equal.
stepUp-monotone : ∀ m run {k} (ps : Vec ℚᵘ k) → Linked _≥_ ps → All (0ℚᵘ ≤_) ps →
                  ∀ (i j : Fin k) → lookup ps i ≤ lookup ps j →
                  lookup (stepUp m run ps) i ≤ lookup (stepUp m run ps) j
stepUp-monotone m run ps lnk all i j pᵢ≤pⱼ with ℕₚ.≤-total (toℕ i) (toℕ j)
... | inj₁ i≤j = ≤-reflexive (stepUp-tie-from m run ps lnk all i j i≤j pᵢ≤pⱼ)
... | inj₂ j≤i = desc-lookup (stepUp-descending m run ps) j i j≤i

------------------------------------------------------------------------
-- `bh_adjust` on the sorted sequence

-- Every adjusted value is ≤ 1.
adjust-≤1 : ∀ {m} (ps : Vec ℚᵘ m) (i : Fin m) → lookup (adjust ps) i ≤ 1ℚᵘ
adjust-≤1 {m} = stepUp-≤run m 1ℚᵘ

-- The adjusted sequence is descending, as the raw one was.
adjust-descending : ∀ {m} (ps : Vec ℚᵘ m) → Linked _≥_ (adjust ps)
adjust-descending {m} = stepUp-descending m 1ℚᵘ

-- q_i ≥ p_i for a descending sequence of p-values in [0, 1].
adjust-≥ : ∀ {m} (ps : Vec ℚᵘ m) → Linked _≥_ ps →
           All (0ℚᵘ ≤_) ps → All (_≤ 1ℚᵘ) ps →
           ∀ (i : Fin m) → lookup ps i ≤ lookup (adjust ps) i
adjust-≥ []       _   _   _         ()
adjust-≥ (p ∷ ps) lnk all (p≤1 ∷ _) i = stepUp-≥ _ 1ℚᵘ p ps ℕₚ.≤-refl p≤1 lnk all i

-- p_i ≤ p_j ⇒ q_i ≤ q_j for a descending sequence of non-negative p-values.
adjust-monotone : ∀ {m} (ps : Vec ℚᵘ m) → Linked _≥_ ps → All (0ℚᵘ ≤_) ps →
                  ∀ (i j : Fin m) → lookup ps i ≤ lookup ps j →
                  lookup (adjust ps) i ≤ lookup (adjust ps) j
adjust-monotone {m} = stepUp-monotone m 1ℚᵘ

------------------------------------------------------------------------
-- `bh_adjust` at the ORIGINAL indices
--
-- `ord = sortperm(q; rev = true)` is a permutation π: sorted position `r`
-- holds the raw p-value at original index `π ⟨$⟩ʳ r`.  `adj[k] = running`
-- writes sorted position `r`'s value back to `ord[r]`, so the adjusted value
-- at original index `k` is the sorted-position value at `π ⟨$⟩ˡ k`.  That π
-- really sorts `p` descending is the named assumption `sorted-desc`; the
-- Agda does not re-implement `sortperm`.

module Unsorted {m} (p : Vec ℚᵘ m) (π : Permutation m m) where

  sorted : Vec ℚᵘ m
  sorted = tabulate (λ r → lookup p (π ⟨$⟩ʳ r))

  -- The adjusted p-value at original index `k`.
  bh : Fin m → ℚᵘ
  bh k = lookup (adjust sorted) (π ⟨$⟩ˡ k)

  -- Reading the sorted sequence at `k`'s sorted position gives `p` at `k`.
  sorted-lookup : ∀ k → lookup sorted (π ⟨$⟩ˡ k) ≡ lookup p k
  sorted-lookup k =
    trans (lookup∘tabulate (λ r → lookup p (π ⟨$⟩ʳ r)) (π ⟨$⟩ˡ k))
          (cong (lookup p) (inverseʳ π))

  -- q_k ≤ 1 at every original index.
  bh-≤1 : ∀ k → bh k ≤ 1ℚᵘ
  bh-≤1 k = adjust-≤1 sorted (π ⟨$⟩ˡ k)

  -- P5, second half: q_k ≥ p_k, given π sorts `p` descending and every
  -- p-value lies in [0, 1].
  bh-≥ : Linked _≥_ sorted →
         (∀ k → 0ℚᵘ ≤ lookup p k) → (∀ k → lookup p k ≤ 1ℚᵘ) →
         ∀ k → lookup p k ≤ bh k
  bh-≥ sorted-desc 0≤p p≤1 k =
    subst (_≤ bh k) (sorted-lookup k)
      (adjust-≥ sorted sorted-desc
        (tabulate⁺ (λ r → 0≤p (π ⟨$⟩ʳ r)))
        (tabulate⁺ (λ r → p≤1 (π ⟨$⟩ʳ r)))
        (π ⟨$⟩ˡ k))

  -- P5, first half: p_i ≤ p_j ⇒ q_i ≤ q_j at any two original indices, given
  -- π sorts `p` descending and every p-value is non-negative.
  bh-monotone : Linked _≥_ sorted → (∀ k → 0ℚᵘ ≤ lookup p k) →
                ∀ i j → lookup p i ≤ lookup p j → bh i ≤ bh j
  bh-monotone sorted-desc 0≤p i j pᵢ≤pⱼ =
    adjust-monotone sorted sorted-desc
      (tabulate⁺ (λ r → 0≤p (π ⟨$⟩ʳ r)))
      (π ⟨$⟩ˡ i) (π ⟨$⟩ˡ j)
      (subst₂ _≤_ (sym (sorted-lookup i)) (sym (sorted-lookup j)) pᵢ≤pⱼ)
