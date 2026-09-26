# Baseline versions and provenance

| Method | Upstream revision | License | Local integration |
|---|---|---|---|
| BBS | v0.2.0, `27c39ec5b631498b3c55073038fa6a59821f2974` | MIT | Unmodified upstream source |
| MUSCLE | v5.3, Bioconda build `h9948957_3` | GNU GPL v3; see `muscle/LICENSE` | External installation plus `run_muscle_consensus.py` |
| CPL | `167a68261719ea458e79c40fcc485ae768c1a64a` | MIT | Upstream source plus `main_cpl32.cpp` |
| ITR | `c50dec739bd2c7f18ac7678d921eecee02e86c6a` | Upstream says `License TBA` | Upstream fetched on demand; local OpenMP adapter |
| TrellisBMA | `77cb3d3965548e0ff63a37ff5de13f96fa15f5a8` | MIT | Required upstream modules plus cluster runner |

## MUSCLE installation provenance

The formal experiments used Bioconda `muscle=5.3`, build `h9948957_3`.
Its executable SHA-256 was:

    1354ce0adce311313ec14574cce64158937ea48f339612ca765802e61626ecfe

An official Linux x86-64 MUSCLE 5.3 binary with build identifier `d9725ac`
was also tested during migration. Its SHA-256 was:

    318abeb951d786a3e2532714cc81ad3b3d8f79a2b517dc31316eeb5b694db2bc

The official binary and the Bioconda build produced byte-identical results
in the 10-cluster migration test.

MUSCLE is installed separately and its executable is not committed to this
repository.

## Verified migration toolchain

- GCC/G++ 9.4.0
- Cargo 1.97.1
- Rust 1.97.1
- TrellisBMA: Python 3.10.20
- MUSCLE wrapper: Python 3.11.15
- MUSCLE: 5.3

## Smoke-test status

All five methods completed the same 10-cluster input with exit status zero
and produced one output sequence per cluster. These integration-test accuracy
values are not paper benchmark results.
