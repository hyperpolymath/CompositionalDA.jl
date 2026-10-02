#!/usr/bin/env bash
# SPDX-License-Identifier: MPL-2.0
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# Bootstrap and run the Agda proof gate for CompositionalDA.jl.
#
# This script exists because the proofs must be checkable by anyone, on a clean
# machine, with one command — and because a proof gate that silently passes when
# the prover is absent is worse than no gate at all.  Every failure mode below
# exits non-zero.  There is no `|| true`, no `command -v agda || exit 0`, and no
# path in which "prover not installed" is reported as success.
#
# Usage:
#   proofs/bootstrap.sh              bootstrap if needed, then check
#   proofs/bootstrap.sh --bootstrap  vendor the standard library only
#   proofs/bootstrap.sh --check      check only (fail if not bootstrapped)
#
# Environment:
#   AGDA_BIN      path to the agda binary to use (default: `agda` on PATH)
#   PROOFS_VENDOR directory for the vendored stdlib / Agda source
#                 (default: proofs/.vendor, which is git-ignored)
#
# This script does NOT install Agda.  The estate bans Python, so the PyPI wheel
# the MetaManifold version of this script fell back to is gone; an absent or
# wrong-version Agda is a loud failure with install instructions (see
# `check_agda_version`).
#
# Pinned versions — these are the versions the proofs were developed against.
# Agda 2.7.0.1 with stdlib 3.0-dev.  Changing either is a real change to the
# gate and must be reviewed, not a maintenance chore to be done silently.

set -euo pipefail

AGDA_VERSION="2.7.0.1"
AGDA_RELEASE_ASSET="Agda-v${AGDA_VERSION}-linux.tar.xz"
AGDA_RELEASE_URL="https://github.com/agda/agda/releases/tag/v${AGDA_VERSION}"

# agda-stdlib revision.  Pinned by SHA, not by tag, and the reason matters:
# these proofs use stdlib APIs (`Data.Integer.Properties._≡?_`, the `Dec` record
# shape, `NonNegative` *instances* for `*-monoʳ-≤-nonNeg`,
# `Algebra.Properties.CommutativeMonoid.Sum.sum-permute`) that exist only on the
# development line towards 3.0 — the branch self-identifies as
# `standard-library-3.0` but no `v3.0` tag has been cut.  They do NOT compile
# against the latest release (v2.4) or against v2.1.  See
# proofs/residue/toolchain.residue.
STDLIB_VERSION="2ffa8b7d4e8e818717ad643d184f055a4d1b0447"

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PROOFS_DIR="$REPO_ROOT/proofs"
AGDA_DIR_SRC="$PROOFS_DIR/agda"
VENDOR="${PROOFS_VENDOR:-$PROOFS_DIR/.vendor}"
LIB_FILE="$AGDA_DIR_SRC/compositionalda.agda-lib"
ENTRY="CompositionalDA/All.agda"

# Print a progress line on stdout, prefixed so it is distinguishable from
# Agda's own output.
log() { printf 'proofs: %s\n' "$*"; }

# Print a fatal error on stderr and exit non-zero.  Every failure in this
# script goes through here; there is no soft-fail path.
die() { printf 'proofs: FATAL: %s\n' "$*" >&2; exit 1; }

[[ -f "$LIB_FILE" ]] || die "missing $LIB_FILE"

##############################################################################
# Locate Agda and insist on the pinned version
##############################################################################

# Set AGDA to the binary named by AGDA_BIN, else `agda` on PATH.  Returns 1 if
# neither exists; AGDA_BIN naming a non-executable is fatal rather than skipped.
resolve_agda() {
  if [[ -n "${AGDA_BIN:-}" ]]; then
    [[ -x "$AGDA_BIN" ]] || die "AGDA_BIN=$AGDA_BIN is not executable"
    AGDA="$AGDA_BIN"
    return 0
  fi
  if command -v agda >/dev/null 2>&1; then
    AGDA="$(command -v agda)"
    return 0
  fi
  return 1
}

# The install instructions printed whenever Agda is absent or the wrong version.
# Written once so the two failure paths cannot drift apart.
install_instructions() {
  cat >&2 <<EOF
proofs: the proof gate needs exactly 'Agda version $AGDA_VERSION'.
proofs: install the official release binary (no package manager, no Python):
proofs:   1. download $AGDA_RELEASE_ASSET from $AGDA_RELEASE_URL
proofs:   2. tar -xJf $AGDA_RELEASE_ASSET -C <dir>
proofs:        this unpacks <dir>/Agda-v$AGDA_VERSION/{bin/agda,data/}
proofs:   3. the release binary hard-codes a CI build path for its data
proofs:      directory, so point it at the bundled one:
proofs:        export Agda_datadir=<dir>/Agda-v$AGDA_VERSION/data
proofs:   4. run the gate with that binary:
proofs:        AGDA_BIN=<dir>/Agda-v$AGDA_VERSION/bin/agda proofs/bootstrap.sh
proofs: (or put a wrapper that sets Agda_datadir on PATH as 'agda').
EOF
}

# Require `agda --version` to print exactly the pinned version.  Anything else —
# a newer Agda, an older one, a binary that cannot run — is fatal.  A proof
# checked by a different prover version than the one it was developed against
# is not the same evidence, and this gate does not pretend otherwise.
check_agda_version() {
  local reported
  reported="$("$AGDA" --version 2>&1 | head -n 1 || true)"
  if [[ "$reported" != "Agda version $AGDA_VERSION" ]]; then
    install_instructions
    die "need 'Agda version $AGDA_VERSION'; '$AGDA --version' printed '${reported:-<nothing>}'"
  fi
}

##############################################################################
# Vendor the standard library and the Agda primitive libraries
##############################################################################

# Fetch agda-stdlib at exactly STDLIB_VERSION into $VENDOR/agda-stdlib and
# normalise its library file to `name: standard-library`, which is what the
# `depend:` clause of compositionalda.agda-lib refers to.  An existing checkout
# is verified against the pin rather than trusted.
vendor_stdlib() {
  if [[ -f "$VENDOR/agda-stdlib/standard-library.agda-lib" ]]; then
    local have
    have="$(git -C "$VENDOR/agda-stdlib" rev-parse HEAD 2>/dev/null || true)"
    [[ "$have" == "$STDLIB_VERSION" ]] \
      || die "$VENDOR/agda-stdlib is at '${have:-<not a git checkout>}', not the pinned $STDLIB_VERSION; remove it and re-run"
    log "stdlib already vendored at $STDLIB_VERSION"
    return 0
  fi
  log "cloning agda-stdlib at $STDLIB_VERSION"
  mkdir -p "$VENDOR"
  rm -rf "$VENDOR/agda-stdlib"
  # A tag would allow `--branch`; a SHA does not, so fetch the exact revision.
  git init -q "$VENDOR/agda-stdlib"
  git -C "$VENDOR/agda-stdlib" remote add origin https://github.com/agda/agda-stdlib
  git -C "$VENDOR/agda-stdlib" fetch -q --depth 1 origin "$STDLIB_VERSION" \
    && git -C "$VENDOR/agda-stdlib" checkout -q FETCH_HEAD \
    || die "could not fetch agda-stdlib $STDLIB_VERSION"
  [[ "$(git -C "$VENDOR/agda-stdlib" rev-parse HEAD)" == "$STDLIB_VERSION" ]] \
    || die "fetched agda-stdlib is not at $STDLIB_VERSION"
  # The upstream `.agda-lib` file is named after the tag (`agda-stdlib.agda-lib`
  # on some tags, `standard-library-2.1.agda-lib` on others, plain
  # `standard-library.agda-lib` at this SHA) and its `name:` does not always read
  # `standard-library` (here it reads `standard-library-3.0`), which is what
  # `defaults` and every `depend:` clause refer to.  Normalise both, whatever the
  # revision ships.
  local lib
  lib="$(find "$VENDOR/agda-stdlib" -maxdepth 1 -name '*.agda-lib' | head -1)"
  [[ -n "$lib" ]] || die "the agda-stdlib checkout at $STDLIB_VERSION has no .agda-lib file"
  sed -i 's/^name: .*/name: standard-library/' "$lib"
  if [[ "$(basename "$lib")" != "standard-library.agda-lib" ]]; then
    cp "$lib" "$VENDOR/agda-stdlib/standard-library.agda-lib"
  fi
  log "stdlib library file: $lib"
}

# Make sure `Agda.Primitive` and friends are available.  The release tarball
# ships them under its data directory (`data/lib/prim`), so with the pinned
# binary this is a no-op.  A binary built without them is fatal: the archive
# version of this script copied them into the data directory, but that
# directory belongs to the toolchain install, not to this repo, and a gate
# must not write into it.
vendor_prims() {
  local datadir
  datadir="$("$AGDA" --print-agda-dir)"
  if [[ -d "$datadir/lib/prim/Agda" ]]; then
    log "primitive libraries present at $datadir/lib/prim"
    return 0
  fi
  install_instructions
  die "$datadir/lib/prim has no Agda primitive libraries; use the release tarball with Agda_datadir pointed at its bundled data/"
}

# The in-repo libraries file that check() passes to Agda via --library-file.
# Only this file is written.  The MetaManifold version of this script also
# copied it into `$(agda --print-agda-dir)/lib` and `$XDG_CONFIG_HOME/agda`;
# the first writes into the toolchain install (which this repo must not touch)
# and neither is needed, because every Agda invocation in this tree names the
# libraries file explicitly.  Interactive users who want `agda` to find the
# libraries without `--library-file` can copy it to `~/.config/agda/libraries`
# themselves.
AGDA_LIBRARIES_FILE="$PROOFS_DIR/.agda-libraries"

# Write proofs/.agda-libraries listing the project library and the vendored
# stdlib, both as absolute paths.
write_agda_config() {
  {
    printf '%s\n' "$LIB_FILE"
    printf '%s\n' "$VENDOR/agda-stdlib/standard-library.agda-lib"
  } > "$AGDA_LIBRARIES_FILE"
  log "wrote $AGDA_LIBRARIES_FILE"
}

##############################################################################
# The gate itself
##############################################################################

# Type-check the whole tree through its single entry module, then run the
# axiom audit.  Both must succeed; the audit is what proves the entry module
# actually reaches every file.
check() {
  resolve_agda || { install_instructions; die "agda is not on PATH and no AGDA_BIN was given"; }
  check_agda_version
  log "using $("$AGDA" --version) at $AGDA"
  cd "$AGDA_DIR_SRC"
  # --safe: no postulates, no foreign code, no `--type-in-type`.  A proof that
  # needs an escape hatch is not a proof.
  # --without-K: no proof-irrelevance-by-fiat for equality.
  log "type-checking CompositionalDA.All (this reaches every module in the tree)"
  [[ -f "$AGDA_LIBRARIES_FILE" ]] \
    || die "$AGDA_LIBRARIES_FILE is missing; run proofs/bootstrap.sh --bootstrap first"
  "$AGDA" --library-file="$AGDA_LIBRARIES_FILE" --safe --without-K "$ENTRY"
  log "axiom audit"
  "$PROOFS_DIR/tests/axiom-audit.sh"
  log "OK"
}

# The bootstrap half: find the prover, vendor the stdlib, write the libraries
# file.  Shared by `--bootstrap` and the no-argument form.
bootstrap() {
  resolve_agda || { install_instructions; die "agda is not on PATH and no AGDA_BIN was given"; }
  check_agda_version
  vendor_stdlib
  vendor_prims
  write_agda_config
}

case "${1:-}" in
  --bootstrap)
    bootstrap
    ;;
  --check)
    check
    ;;
  "")
    bootstrap
    check
    ;;
  *)
    die "unknown argument: $1 (expected --bootstrap, --check, or nothing)"
    ;;
esac
