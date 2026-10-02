<!-- SPDX-License-Identifier: CC-BY-SA-4.0 -->
<!-- SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) -->

# Oracle datasets

The R-reference oracle (`test/oracle/runoracle.jl`) compares each native method
with its Bioconductor implementation on three two-group count tables.

| Name | Kind | Samples | Features | Groups | Known answer |
|---|---|---|---|---|---|
| `mock` | synthetic, generated in `test/oracle/datasets.jl` | 20 | 40 | a 10 / b 10 | yes: features 1–3 up 8×, 4–5 down 4× |
| `gut` | stool 16S, Parkinson's disease vs healthy controls | 108 | 202 | H 55 / PAR 53 | no |
| `soil` | blueberry soil vs rhizosphere 16S | 63 | 300 | Soil 39 / Rhizosphere 24 | no |

`mock` is not a MetaManifold fixture; that repository holds no count-table
fixtures. Real tables have no ground truth, so for `gut` and `soil` the oracle
checks agreement with R only.

## Source and licence of `gut` and `soil`

Both come from the deposit accompanying Nearing et al. (2022), which compared
differential-abundance methods on 38 two-group 16S datasets:

- Nearing J (2021). *16S_rRNA_Microbiome_Datasets*. figshare.
  doi:[10.6084/m9.figshare.14531724.v1](https://doi.org/10.6084/m9.figshare.14531724.v1).
  Licence **CC BY 4.0** (as stated by the figshare record).
- Nearing JT et al. (2022). Microbiome differential abundance methods produce
  different results across 38 datasets. *Nat Commun* 13:342.
  doi:10.1038/s41467-022-28034-z.

File: `DA_COMPARE_DATA_21_03_08.tar.gz`, 253,151,868 bytes, md5
`3bf04686142ba000b3b6f7825cc1a606` (as published by figshare), sha256
`53337642e26735fc138ed17089f6a467135719413bafbd1f9e5bff048c164b84`.
Members used: `Hackathon/Studies/par_scheperjans/` and
`Hackathon/Studies/Blueberry/`, the raw `*_ASVs_table.tsv` (not the rarefied
`*_rare.tsv`) and the metadata TSV (`comparison` column).

**Changes made** (CC BY 4.0 §3(a)(1)(B)): samples intersected with the
metadata; features kept if present in ≥ 25% of samples, then at most the 300
most prevalent; transposed to the package's samples × features convention at
load time. `test/oracle/make_fixtures.jl` reproduces the fixtures exactly from
the checksummed tarball. Each fixture carries a REUSE `.license` sidecar.

The deposit's CC BY 4.0 covers these processed tables. The licences of the
original studies behind each table were not separately verified, and the
original-study citations are not recorded here because they were not confirmed
from the deposit itself.

The fixture data are CC BY 4.0, not MPL-2.0; they are test data only and are
not part of the package's loaded code.
