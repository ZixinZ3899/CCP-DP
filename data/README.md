# CCP-DP datasets

Large sequencing datasets are not committed to this repository. This directory
records the exact benchmark files, their provenance, and reproducible preparation
steps where the complete transformation chain is known.

## Dataset inventory

| ID | Centers | Length | Role | Reproduction status |
|---|---:|---:|---|---|
| `srinivas_110` | 9,984 | 110 bp | Primary benchmark | Fully reproducible |
| `chandak_108` | 1,466 | 108 bp | Primary benchmark | Downloadable processed benchmark |
| `bar_lev_140` | 10,000 | 140 bp | Primary benchmark | Downloadable subsampled benchmark |
| `chandak_150` | 11,710 | 150 bp | Additional real-data benchmark | Source public; local benchmark conversion partly documented |
| `erlich_152` | 72,000 | 152 bp | Additional real-data benchmark | Raw reads public; benchmark clustering pipeline unavailable |
| `grass_117` | 4,982 | 117 bp | Archived additional benchmark | Local processed benchmark; conversion pipeline unavailable |

Machine-readable paths, checksums, citations, and source URLs are in
`datasets.json`. Expected SHA-256 values are in `SHA256SUMS`.

## Quick setup: primary benchmarks

From the repository root:

```bash
./data/download_datasets.sh --core
python3 data/prepare_datasets.py validate --config data/datasets.json
```

This creates:

```text
data/processed/
├── srinivas_110/{Clusters.txt,Centers.txt}
├── chandak_108/{Clusters.txt,Centers.txt}
└── bar_lev_140/{Clusters.txt,Centers.txt}
```

The downloader does not place large data files under version control.

## Primary datasets

### Srinivas-110

The original Microsoft CNR repository contains 10,000 centers of length 110
and a cluster file that includes empty clusters. CCP-DP removes the 16 empty
clusters and their aligned centers, preserving the order of all remaining
clusters. The resulting benchmark contains 9,984 cluster-center pairs.

- Repository: <https://github.com/microsoft/clustered-nanopore-reads-dataset>
- Pinned commit: `6938f44796185902a08381943c2895782886c5c3`
- Original files: `Clusters.txt`, `Centers.txt`
- Preparation: `prepare_datasets.py filter-empty`

The upstream repository notes that the generated centers exhibit unintended
long-range dependencies. This upstream caveat must be retained when describing
the dataset.

### Chandak-108

CCP-DP uses the processed trace-reconstruction benchmark deposited in Zenodo
record 16959565:

- `oligo0_UnderlyingClusters.txt` -> `Clusters.txt`
- `oligo0refs.txt` -> `Centers.txt`

The Zenodo record states that this version uses perfect clustering: sequenced
reads are assigned to their closest ground-truth sequence. It is therefore a
processed benchmark derived from the Chandak et al. experiment, not an
unmodified raw-data download.

- Benchmark record: <https://zenodo.org/records/16959565>
- Original data repository: <https://github.com/shubhamchandak94/nanopore_dna_storage_data>

### Bar-Lev-140

CCP-DP uses the 10,000-cluster subsampled benchmark from Zenodo record 16959565:

- `BinnedNanoporeTwoFlowcells_clusters_subsampled.txt` -> `Clusters.txt`
- `BinnedNanoporeTwoFlowcells_centers_subsampled.txt` -> `Centers.txt`

This is a post-processed and subsampled derivative of the Bar-Lev et al.
dataset, not the complete original dataset. The Zenodo record identifies the
full source as record 13896773.

## Additional datasets

### Chandak-150

The original code and data are available from:

- <https://github.com/shubhamchandak94/LDPC_DNA_storage>
- <https://github.com/shubhamchandak94/LDPC_DNA_storage_data>

The CCP-DP benchmark uses `Chandak_ref28_Centers.txt` and a corresponding
cluster file. The local `_31eq` derivative changes only separator lines to 31
equals signs; it does not alter read sequences or cluster order:

```bash
python3 data/prepare_datasets.py normalize-separators \
  --input PATH/TO/Chandak_ref28_Clusters.txt \
  --output data/processed/chandak_150/Clusters.txt
```

The earlier steps that produced the `ref28` benchmark pair were not retained,
so this dataset is not downloaded or claimed as end-to-end reproducible by the
default script.

### Erlich-152

The DNA Fountain repository documents raw sequencing accessions
`PRJEB19305` and `PRJEB19307` and the original read-processing procedure:

- <https://github.com/TeamErlich/dna-fountain>
- <https://www.ebi.ac.uk/ena/browser/view/PRJEB19305>
- <https://www.ebi.ac.uk/ena/browser/view/PRJEB19307>

The exact transformation from those raw reads to the local 72,000
cluster-center benchmark was not found. Consequently the repository records
checksums for the benchmark but does not provide a misleading automatic
reconstruction recipe.

### Grass-117

The Grass et al. experiment used 158-nt oligonucleotides containing a 117-nt
data payload plus sequencing adapters. The local benchmark contains 4,982
117-bp centers, whereas the paper reports 4,991 synthesized segments. No script
explaining the removal of nine entries was found. This pair is therefore
treated as a preserved processed benchmark, not as an automatically generated
copy of the publication data.

## Installing preserved manual benchmark files

If the six additional files from the audited legacy workspace are available,
place them as follows:

```text
data/processed/chandak_150/Clusters.txt
data/processed/chandak_150/Centers.txt
data/processed/erlich_152/Clusters.txt
data/processed/erlich_152/Centers.txt
data/processed/grass_117/Clusters.txt
data/processed/grass_117/Centers.txt
```

Then validate all present files:

```bash
python3 data/prepare_datasets.py validate --config data/datasets.json
```

Missing optional datasets are reported as `SKIP`; a present file with the wrong
checksum is an error.

## Data reporting rules

1. Cite the original experimental paper for every dataset used.
2. Cite Zenodo record 16959565 when using its processed Chandak or Bar-Lev files.
3. Describe Chandak-108 as perfectly clustered and Bar-Lev-140 as subsampled.
4. Describe Srinivas-110 as the non-empty 9,984-cluster derivative.
5. Do not describe Erlich-152 or Grass-117 as end-to-end regenerated datasets.
6. Downsampled coverage experiments must save their cluster indices, random
   seeds, and generated center files together with their results.

Full bibliographic entries are provided in `CITATIONS.bib`.
