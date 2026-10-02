#!/usr/bin/env bash
# SPDX-License-Identifier: MPL-2.0
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# The whole proof gate, in the order CI runs it:
#   1. axiom audit (--safe on every module, no postulates, FFI, unsound flags,
#      trusted primitives or holes; every module reachable from All.agda)
#   2. type-check CompositionalDA/All.agda; any Agda warning is a failure
#   3. gate self-test: the three reject/ controls must fail for their stated
#      reason, and each audit check must fire on a planted defect
#
# Toolchain: Debian 13's agda (2.6.4.3) and agda-stdlib (2.1), resolved by
# proofs/lib.sh.  Nothing is downloaded or vendored.
#
# Usage: proofs/check.sh

set -euo pipefail
# shellcheck source=proofs/lib.sh
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
resolve_agda
resolve_stdlib
printf 'check: %s, %s\n' "$("$AGDA" --version)" "$STDLIB_PINNED_NAME"

"$PROOFS_DIR/tests/axiom-audit.sh"

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
write_libraries "$WORK/libraries" "$AGDA_DIR"
# --safe and --without-K come from compositionalda.agda-lib, so a module cannot
# opt out of them; they are repeated here so the command line says what it does.
( cd "$AGDA_DIR" && "$AGDA" --library-file="$WORK/libraries" --safe --without-K "$ENTRY" ) 2>&1 | tee "$WORK/typecheck.log"
if grep -qiE '^warning|^[^ ].*: *warning' "$WORK/typecheck.log"; then
  die "Agda emitted warnings"
fi
printf 'check: %s type-checks\n' "$ENTRY"

"$PROOFS_DIR/tests/gate-selftest.sh"
printf 'check: OK\n'
