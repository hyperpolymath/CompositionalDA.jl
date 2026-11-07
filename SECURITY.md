<!-- SPDX-License-Identifier: CC-BY-SA-4.0 -->
<!-- SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) -->

# Security and numerical integrity

This is a development package; no released version currently carries a
security support commitment. Report exploitable defects through GitHub's
private vulnerability-reporting facility if enabled. Otherwise contact the
maintainer at the address in Project.toml before publishing exploit details.
Ordinary numerical defects can be reported as issues with a minimal synthetic
reproducer, package/tool versions and expected versus observed results.
Do not attach confidential or identifiable biological data.

The native API is a local scientific library, not a sandbox for untrusted jobs.
Callers must bound matrix dimensions, count bit widths, Monte-Carlo work,
concurrency and memory at their service boundary. Input validation and explicit
nonconvergence errors do not create a hard CPU or memory limit. R reference
execution loads installed R packages and must run in a separately controlled
environment when inputs or dependencies are untrusted.

Aqua checks Julia package quality. It is not a vulnerability-database scanner.
The launch ledger must separately record dependency advisories, secret scanning,
workflow permissions/action pins, release provenance and unresolved findings.
No successful comprehensive security audit is claimed here. Statistics that
fail calibration are safety/integrity concerns even without an attacker.
