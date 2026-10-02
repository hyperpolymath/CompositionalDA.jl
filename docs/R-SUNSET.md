<!-- SPDX-License-Identifier: CC-BY-SA-4.0 -->
<!-- SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) -->

# Sunset criterion for the R reference backend

> **Status:** the extension described here is planned and not yet in the tree:
> there is no `ext/` directory and no `test/oracle_*.jl` yet.

The R backend (`ext/CompositionalDARCallExt.jl`, loaded only when RCall is
present) is a bridge, not a product. It never owns R state. Callers pass an
`evaluator` function, so a host application keeps its own R session and lock.

**Criterion.** For each method, once the native implementation agrees with the R
reference within the tolerances stated in `test/oracle_*.jl` on all three
reference datasets (mock, gut, soil) for **two consecutive releases**, that
method's R path is marked deprecated in the next release and receives no
further maintenance. When no method still needs it, the extension is deleted.
That removal touches only `ext/` and the `[weakdeps]`/`[extensions]` sections of
`Project.toml`.

**What stays after sunset.** The recorded oracle outputs (fixtures under
`test/fixtures/`) remain as frozen regression data. Agreement with them is then
tested without R being installed.
