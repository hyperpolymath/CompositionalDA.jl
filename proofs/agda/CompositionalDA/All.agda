-- SPDX-License-Identifier: MPL-2.0
-- SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
--
-- The gate entry point: type-checking this module type-checks every module in
-- `proofs/agda/CompositionalDA`.  `proofs/tests/axiom-audit.sh` verifies that
-- every module in the tree is reachable from here, so a module that is not
-- imported below is an audit failure, not a silent omission.

{-# OPTIONS --without-K --safe #-}

module CompositionalDA.All where

import CompositionalDA.Prelude
import CompositionalDA.LogHom
import CompositionalDA.CLR
import CompositionalDA.CLR.Instance
import CompositionalDA.BenjaminiHochberg
import CompositionalDA.BenjaminiHochberg.StepUp
