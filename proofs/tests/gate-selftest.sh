#!/usr/bin/env bash
# SPDX-License-Identifier: MPL-2.0
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# Show that the proof gate rejects what it claims to reject.
#
# A gate that has never been seen to fail is not evidence of anything.  This
# script plants known defects and requires the gate to catch every one:
#
#   A. Each negative control in proofs/agda/reject/ must FAIL to type-check,
#      and Agda's output must match the control's own `-- EXPECT: <regex>`
#      line.  A control that fails for some other reason (a missing import, a
#      typo) does not count: it would keep "passing" after the defect it
#      guards against had been reintroduced.
#   B. Each of the seven checks in axiom-audit.sh must flag a defect planted
#      in a scratch copy of the proof tree.  The unmodified copy must pass
#      first, so a planted failure cannot be an artefact of the copy.
#
# Prerequisite: proofs/bootstrap.sh has been run (it writes the libraries file
# and vendors the standard library).  Every failure exits non-zero.

set -euo pipefail

PROOFS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
AGDA_DIR="$PROOFS_DIR/agda"
LIBRARIES_FILE="$PROOFS_DIR/.agda-libraries"
AGDA="${AGDA_BIN:-agda}"

failures=0

# Record one finding on stderr and keep going, so a run reports every gap.
fail() { printf 'gate-selftest: FAIL: %s\n' "$*" >&2; failures=$((failures + 1)); }

# Print a progress line on stdout.
note() { printf 'gate-selftest: %s\n' "$*"; }

[[ -f "$LIBRARIES_FILE" ]] \
  || { printf 'gate-selftest: FATAL: %s missing; run proofs/bootstrap.sh first\n' "$LIBRARIES_FILE" >&2; exit 1; }
command -v "$AGDA" >/dev/null 2>&1 || [[ -x "$AGDA" ]] \
  || { printf 'gate-selftest: FATAL: agda not found (set AGDA_BIN)\n' >&2; exit 1; }

SCRATCH="$(mktemp -d)"
trap 'rm -rf "$SCRATCH"' EXIT

##############################################################################
# A. Negative controls must fail, for the stated reason
##############################################################################

# Type-check one reject/ control and require a non-zero exit whose output
# matches the regex on the control's `-- EXPECT:` line.
check_reject() {
  local file="$1" rel expect out rc
  rel="${file#"$AGDA_DIR"/}"
  expect="$(sed -n 's/^-- EXPECT: //p' "$file" | head -1)"
  if [[ -z "$expect" ]]; then
    fail "$rel has no '-- EXPECT:' line, so its failure reason cannot be checked"
    return
  fi
  set +e
  out="$(cd "$AGDA_DIR" && "$AGDA" --library-file="$LIBRARIES_FILE" --safe --without-K "$rel" 2>&1)"
  rc=$?
  set -e
  if [[ $rc -eq 0 ]]; then
    fail "$rel type-checked; the gate accepts a defect it must reject"
  elif ! grep -Eq -- "$expect" <<<"$out"; then
    fail "$rel failed, but not for the expected reason /$expect/: $(tail -3 <<<"$out" | tr '\n' ' ')"
  else
    note "rejected as expected: $rel"
  fi
}

mapfile -t REJECTS < <(find "$AGDA_DIR/reject" -name '*.agda' -type f | sort)
[[ ${#REJECTS[@]} -gt 0 ]] || fail "no negative controls under $AGDA_DIR/reject"
for f in "${REJECTS[@]}"; do check_reject "$f"; done

##############################################################################
# B. Every axiom-audit check must fire on a planted defect
##############################################################################

# Make a fresh scratch copy of the audited tree (agda/ and tests/) in $1.
fresh_copy() {
  rm -rf "$1"
  mkdir -p "$1"
  cp -r "$AGDA_DIR" "$1/agda"
  cp -r "$PROOFS_DIR/tests" "$1/tests"
  rm -rf "$1/agda/_build"
}

# Run the audit inside the scratch copy $1; prints its combined output and
# returns its exit status.
run_audit() { "$1/tests/axiom-audit.sh" 2>&1; }

CLEAN="$SCRATCH/clean"
fresh_copy "$CLEAN"
if ! out="$(run_audit "$CLEAN")"; then
  fail "the unmodified copy fails the audit, so planted results would mean nothing: $out"
else
  note "positive control: the unmodified copy passes the audit"
fi

TARGET="agda/CompositionalDA/Prelude.agda"

# Plant one defect with the shell snippet $2 (run inside the copy), then
# require the audit to exit non-zero with a finding matching the regex $3.
plant() {
  local name="$1" mutate="$2" expect="$3" dir="$SCRATCH/$1" out
  fresh_copy "$dir"
  (cd "$dir" && eval "$mutate")
  if out="$(run_audit "$dir")"; then
    fail "audit check '$name' missed its planted defect"
  elif ! grep -Eq -- "$expect" <<<"$out"; then
    fail "audit check '$name' failed, but without the expected finding /$expect/: $out"
  else
    note "audit caught planted defect: $name"
  fi
}

plant safe      "sed -i 's/ --safe//' $TARGET"                                   'does not declare --safe'
plant postulate "printf '\npostulate\n  bogus : Set\n' >> $TARGET"               'postulate block'
plant ffi       "printf '\n{-# FOREIGN GHC import Data.List #-}\n' >> $TARGET"   'FOREIGN/COMPILE/BUILTIN'
plant unsound   "sed -i 's/--safe/--safe --type-in-type/' $TARGET"               'unsound flag'
plant trust     "printf '\n-- uses trustMe\n' >> $TARGET"                        'trustMe/primTrust'
plant hole      "printf '\nhole = {! !}\n' >> $TARGET"                           'contains a hole'
plant reach     "printf '{-# OPTIONS --without-K --safe #-}\nmodule CompositionalDA.Orphan where\n' > agda/CompositionalDA/Orphan.agda" \
                'not reachable'

if [[ $failures -gt 0 ]]; then
  printf 'gate-selftest: %d gap(s) in the gate\n' "$failures" >&2
  exit 1
fi
note "OK: ${#REJECTS[@]} negative control(s) rejected, 7 audit checks fire"
