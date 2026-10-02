#!/usr/bin/env bash
# SPDX-License-Identifier: MPL-2.0
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# Prove that the proof gate can fail.
#
# A gate that has never been observed to reject anything is not evidence.  This
# script takes a throwaway copy of the proof tree, breaks it in one specific way
# at a time, and requires the gate to reject it.  If any mutation is *accepted*,
# or if a mutation failed to apply at all, this script exits non-zero.
#
# Three kinds of control:
#   1. planted defects in a theorem — the type-checker must reject them;
#   2. planted defects in the tree's hygiene (a module dropped from the entry
#      point, a postulate, `--safe` removed) — the axiom audit must reject them;
#   3. the standing negative controls in `proofs/agda/reject/` — Agda must
#      reject each one with the exact error its `-- EXPECT:` line names, so that
#      a rejection for the wrong reason (a missing import, say) does not pass
#      for the right one.
# And finally the pristine tree must still be accepted, so that a gate which
# rejects everything cannot pass this test either.
#
# It runs in CI.  A gate whose self-test is skipped is a gate nobody can trust,
# so there is no flag to turn this off.
#
# Usage: proofs/tests/gate-selftest.sh
#
# Environment (same interface as proofs/bootstrap.sh):
#   AGDA_BIN      path to the agda binary (default: `agda` on PATH)
#   PROOFS_VENDOR where bootstrap.sh vendored the stdlib (default: proofs/.vendor)

set -uo pipefail

PROOFS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
NAMESPACE="CompositionalDA"
LIB_BASENAME="compositionalda.agda-lib"
VENDOR="${PROOFS_VENDOR:-$PROOFS_DIR/.vendor}"

# Print a fatal error on stderr and exit 1.  Every setup failure goes through
# here; a self-test that cannot run must not look like a self-test that passed.
die() { printf 'gate-selftest: FATAL: %s\n' "$*" >&2; exit 1; }

# Set AGDA from AGDA_BIN or PATH, the same way bootstrap.sh does.  There is no
# vendored-binary fallback: the estate installs Agda from the release tarball,
# and an absent prover is a failure here, not a skip.
resolve_agda() {
  if [[ -n "${AGDA_BIN:-}" ]]; then
    AGDA="$AGDA_BIN"
  elif command -v agda >/dev/null 2>&1; then
    AGDA="$(command -v agda)"
  else
    die "agda not found (run proofs/bootstrap.sh, or set AGDA_BIN)"
  fi
  [[ -x "$AGDA" ]] || die "$AGDA is not executable"
}

# Set STDLIB_LIBS to the library files the gate depends on, read from the
# `proofs/.agda-libraries` that bootstrap.sh wrote (every line except the
# project's own library file).  Falls back to the conventional vendored
# location if that file is absent, and fails if neither exists.
resolve_stdlib() {
  STDLIB_LIBS=()
  local libs_file="$PROOFS_DIR/.agda-libraries" line
  if [[ -f "$libs_file" ]]; then
    while IFS= read -r line; do
      [[ -z "$line" ]] && continue
      [[ "$(basename "$line")" == "$LIB_BASENAME" ]] && continue
      [[ -f "$line" ]] || die "$libs_file names $line, which does not exist"
      STDLIB_LIBS+=("$line")
    done < "$libs_file"
  elif [[ -f "$VENDOR/agda-stdlib/standard-library.agda-lib" ]]; then
    STDLIB_LIBS+=("$VENDOR/agda-stdlib/standard-library.agda-lib")
  fi
  [[ ${#STDLIB_LIBS[@]} -gt 0 ]] \
    || die "standard library not found; run proofs/bootstrap.sh first"
}

resolve_agda
resolve_stdlib

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

# Replace the working copy with a fresh copy of the pristine proof tree and
# the test scripts.  Interface files (`_build/`, `*.agdai`) are removed from
# the copy so that nothing checked against the pristine sources can be reused
# for a mutated one.
reset_tree() {
  rm -rf "$WORK/agda" "$WORK/tests"
  cp -r "$PROOFS_DIR/agda" "$WORK/agda"
  cp -r "$PROOFS_DIR/tests" "$WORK/tests"
  rm -rf "$WORK/agda/_build"
  find "$WORK/agda" -name '*.agdai' -delete
  [[ -f "$WORK/agda/$LIB_BASENAME" ]] || die "missing $PROOFS_DIR/agda/$LIB_BASENAME"
}

reset_tree

# Point Agda at the *copy*: without this the library file would resolve
# `CompositionalDA.All` back to the pristine tree and every mutation below would
# be checked against unmutated sources — a self-test that tests nothing.
{
  printf '%s\n' "$WORK/agda/$LIB_BASENAME"
  printf '%s\n' "${STDLIB_LIBS[@]}"
} > "$WORK/libraries"

# Type-check the copy's entry module with the gate's exact flags.  Output goes
# to $WORK/out.txt; the exit status is the verdict.
run_gate() {
  ( cd "$WORK/agda" && "$AGDA" --library-file="$WORK/libraries" --safe --without-K "$NAMESPACE/All.agda" ) \
    >"$WORK/out.txt" 2>&1
}

# Run the copied axiom audit, which resolves its own tree from its location
# and so audits the copy.
run_audit() { "$WORK/tests/axiom-audit.sh" >"$WORK/out.txt" 2>&1; }

total=0
failed=0

# expect_reject <name> <runner> <file> <mutation-script> [reason-regex]
#
# Restore <file> (relative to proofs/agda) in the copy, apply the mutation
# (a shell snippet receiving the file as $1), and require <runner> to fail.
# Counts a FAIL if the mutation did not change the file (a stale pattern
# would otherwise pass silently as "rejected") or if the runner accepted the
# mutated tree.  With a [reason-regex], the runner must also fail WITH that
# finding: a rejection for some other reason is a FAIL.  Restores the file
# afterwards.
expect_reject() {
  local name="$1" runner="$2" target="$3" mutator="$4" reason="${5:-}"
  local file="$WORK/agda/$target"
  total=$((total + 1))

  cp "$PROOFS_DIR/agda/$target" "$file"
  local before after
  before="$(sha256sum "$file" | awk '{print $1}')"
  bash -c "$mutator" mutator "$file"
  after="$(sha256sum "$file" | awk '{print $1}')"
  if [[ "$before" == "$after" ]]; then
    printf 'gate-selftest: FAIL  %-48s mutation did not apply (stale pattern?)\n' "$name"
    failed=$((failed + 1))
    cp "$PROOFS_DIR/agda/$target" "$file"
    return
  fi

  if $runner; then
    printf 'gate-selftest: FAIL  %-48s gate ACCEPTED a broken proof\n' "$name"
    sed 's/^/                     | /' "$WORK/out.txt" | head -6
    failed=$((failed + 1))
  elif [[ -n "$reason" ]] && ! grep -qE -- "$reason" "$WORK/out.txt"; then
    printf 'gate-selftest: FAIL  %-48s rejected, but not for the expected reason\n' "$name"
    printf '                     | expected /%s/\n' "$reason"
    sed 's/^/                     | /' "$WORK/out.txt" | head -6
    failed=$((failed + 1))
  else
    printf 'gate-selftest: ok    %-48s rejected\n' "$name"
  fi
  cp "$PROOFS_DIR/agda/$target" "$file"
}

# expect_reject_for_reason <file>
#
# Type-check one standing negative control from proofs/agda/reject/ and
# require BOTH that Agda exits 42 (a type error, as opposed to 1 for a
# missing file or library) AND that its output matches the extended regex on
# the control's `-- EXPECT:` line.  A control with no EXPECT line is a FAIL:
# an unexplained rejection is not a control.
expect_reject_for_reason() {
  local target="$1" name file pattern rc
  name="reject/$(basename "${target%.agda}")"
  file="$WORK/agda/$target"
  total=$((total + 1))

  pattern="$(sed -nE 's/^-- EXPECT: (.*)$/\1/p' "$file" | head -1)"
  if [[ -z "$pattern" ]]; then
    printf 'gate-selftest: FAIL  %-48s has no -- EXPECT: line\n' "$name"
    failed=$((failed + 1))
    return
  fi

  ( cd "$WORK/agda" && "$AGDA" --library-file="$WORK/libraries" --safe --without-K "$target" ) \
    >"$WORK/out.txt" 2>&1
  rc=$?
  if [[ $rc -eq 0 ]]; then
    printf 'gate-selftest: FAIL  %-48s gate ACCEPTED a negative control\n' "$name"
    failed=$((failed + 1))
  elif [[ $rc -ne 42 ]]; then
    printf 'gate-selftest: FAIL  %-48s exited %d, not 42 (not a type error)\n' "$name" "$rc"
    sed 's/^/                     | /' "$WORK/out.txt" | head -6
    failed=$((failed + 1))
  elif ! grep -qE -- "$pattern" "$WORK/out.txt"; then
    printf 'gate-selftest: FAIL  %-48s rejected, but not for the expected reason\n' "$name"
    printf '                     | expected /%s/\n' "$pattern"
    sed 's/^/                     | /' "$WORK/out.txt" | grep -v '^ *| *Checking' | head -6
    failed=$((failed + 1))
  else
    printf 'gate-selftest: ok    %-48s rejected for the expected reason\n' "$name"
  fi
}

# --- the type-checker must reject wrong mathematics -------------------------
#
# Each mutator is a shell snippet receiving the target file as $1.  If a snippet
# stops matching (because the proof was rewritten), `expect_reject` reports the
# mutation as unapplied rather than silently passing.

expect_reject "clr: sum-to-zero claimed as sum-to-one" run_gate "$NAMESPACE/CLR.agda" \
  'sed -i "s|clr-sum-zero : ∀ (x : Vector Pos n) → sum (clr x) ≈ 0#|clr-sum-zero : ∀ (x : Vector Pos n) → sum (clr x) ≈ 1#|" "$1"'

expect_reject "clr: centring removed from the transform" run_gate "$NAMESPACE/CLR.agda" \
  'sed -i "s|^clr x = centre (map log x)|clr x = map log x|" "$1"'

expect_reject "bh: family-size monotonicity reversed" run_gate "$NAMESPACE/BenjaminiHochberg.agda" \
  'sed -i "s|  bhScale M j n d ≤ℚ bhScale M′ j n d|  bhScale M′ j n d ≤ℚ bhScale M j n d|" "$1"'

# Mirrors test/mutants.jl `mutant_bh_no_cummin`: scale each rank, never take
# the running minimum.
expect_reject "stepup: running minimum dropped" run_gate "$NAMESPACE/BenjaminiHochberg/StepUp.agda" \
  'sed -i "s|  (run ⊓ scale m k p) ∷ stepUp m (run ⊓ scale m k p) ps|  scale m k p ∷ stepUp m run ps|" "$1"'

expect_reject "stepup: q ≥ p reversed" run_gate "$NAMESPACE/BenjaminiHochberg/StepUp.agda" \
  'sed -i "s|           ∀ (i : Fin m) → lookup ps i ≤ lookup (adjust ps) i|           ∀ (i : Fin m) → lookup (adjust ps) i ≤ lookup ps i|" "$1"'

expect_reject "stepup: monotonicity reversed" run_gate "$NAMESPACE/BenjaminiHochberg/StepUp.agda" \
  'sed -i "s|                  lookup (adjust ps) i ≤ lookup (adjust ps) j|                  lookup (adjust ps) j ≤ lookup (adjust ps) i|" "$1"'

expect_reject "stepup: scale written as rank/m" run_gate "$NAMESPACE/BenjaminiHochberg/StepUp.agda" \
  'sed -i "s|^scale m k p = p \* mkℚᵘ (+ m) k|scale m k p = p * mkℚᵘ (+ suc k) m|" "$1"'

# --- the axiom audit must reject an unchecked or unsound module -------------

expect_reject "audit: module dropped from the gate entry" run_audit "$NAMESPACE/All.agda" \
  'sed -i "/import CompositionalDA.BenjaminiHochberg.StepUp/d" "$1"' 'not reachable'

expect_reject "audit: postulate injected" run_audit "$NAMESPACE/BenjaminiHochberg/StepUp.agda" \
  'printf "postulate cheat : ∀ {A : Set} → A\n" >> "$1"' 'postulate block'

expect_reject "audit: --safe removed" run_audit "$NAMESPACE/CLR.agda" \
  'sed -i "s|--without-K --safe|--without-K|" "$1"' 'does not declare --safe'

expect_reject "audit: FFI pragma injected" run_audit "$NAMESPACE/Prelude.agda" \
  'printf "\n{-# FOREIGN GHC import Data.List #-}\n" >> "$1"' 'FOREIGN/COMPILE/BUILTIN'

expect_reject "audit: unsound flag added" run_audit "$NAMESPACE/Prelude.agda" \
  'sed -i "s|--safe|--safe --type-in-type|" "$1"' 'unsound flag'

expect_reject "audit: trustMe mentioned" run_audit "$NAMESPACE/Prelude.agda" \
  'printf "\n-- uses trustMe\n" >> "$1"' 'trustMe/primTrust'

expect_reject "audit: hole left in a module" run_audit "$NAMESPACE/Prelude.agda" \
  'printf "\nhole = {! !}\n" >> "$1"' 'contains a hole'

# --- the standing negative controls must fail for their stated reason ------

while IFS= read -r control; do
  expect_reject_for_reason "reject/$(basename "$control")"
done < <(find "$PROOFS_DIR/agda/reject" -maxdepth 1 -name '*.agda' -type f | sort)

# --- the gate must still accept the pristine tree --------------------------

total=$((total + 1))
reset_tree
if run_gate && run_audit; then
  printf 'gate-selftest: ok    %-48s accepted\n' "pristine tree"
else
  printf 'gate-selftest: FAIL  %-48s gate REJECTED the real proofs\n' "pristine tree"
  sed 's/^/                     | /' "$WORK/out.txt" | head -12
  failed=$((failed + 1))
fi

printf 'gate-selftest: %d/%d controls behaved correctly\n' "$((total - failed))" "$total"
[[ $failed -eq 0 ]] || exit 1
