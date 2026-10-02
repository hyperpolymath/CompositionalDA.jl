-- SPDX-License-Identifier: MPL-2.0
-- SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
--
-- Negative control: a postulate must be rejected under --safe.
--
-- No OPTIONS pragma here on purpose: the `--safe --without-K` flags come from
-- `compositionalda.agda-lib` and from the gate's command line, and this file
-- checks that those flags actually bite.  `proofs/tests/axiom-audit.sh` also
-- greps for `postulate` independently of Agda.
--
-- EXPECT: [Pp]ostulate

module reject.Postulate where

open import Data.Nat.Base using (ℕ)
open import Relation.Binary.PropositionalEquality using (_≡_)

postulate every-natural-is-zero : (n : ℕ) → n ≡ 0
