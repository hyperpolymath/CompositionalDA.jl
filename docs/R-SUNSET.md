<!-- SPDX-License-Identifier: CC-BY-SA-4.0 -->
<!-- SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) -->

# Sunset criterion for the R reference backend

> **Status (v0.3 development):** the previous 13-check ALDEx2 oracle did not
> cover descriptive-column units or identical-instance q-values. The stronger
> gate must establish a new baseline before a release qualifies for sunset.
> ANCOM-BC now has a reference adapter and conformance gate; their existence
> does not mean the gate has passed. Frozen reference-output fixtures under
> `test/fixtures/` do not yet exist.

The R backend (`ext/CompositionalDARCallExt.jl`, loaded only when RCall is
present) is a bridge, not a product. It never owns R state. Every R call runs inside a
caller-supplied `runner` function, so a host application keeps its own R session
and lock (MetaManifold passes its `with_r_lock`).

**Criterion.** For each implemented method, once the native implementation agrees with the R reference
within the tolerances stated in `test/oracle/runoracle.jl` on all three
reference datasets (mock, gut, soil) for **two consecutive releases**, that
method's R path is marked deprecated in the next release and receives no
further maintenance. When no method still needs it, the extension is deleted.
That removal touches only `ext/` and the `[weakdeps]`/`[extensions]` sections of
`Project.toml`.

**What stays after sunset.** The recorded oracle outputs (fixtures under
`test/fixtures/`) remain as frozen regression data. Agreement with them is then
tested without R being installed.
