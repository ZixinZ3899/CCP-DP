# CCP-DP reproducibility record

This document records the clean-clone reproduction of the main CCP-DP
implementation performed on 26 September 2026.

## Reproduction scope

The validation used a fresh clone of the public Git tree rather than the
original development working directory.

The following operations were performed:

1. clone the public repository history;
2. verify the prepared dataset checksums;
3. compile `src/main.cpp` from scratch;
4. reconstruct every cluster in five complete datasets;
5. verify output counts and target sequence lengths;
6. evaluate every reconstruction against the corresponding centers;
7. compare the primary results with the preserved formal results;
8. record output hashes and the software and hardware environment.

This was a full-dataset reproduction, not a limited smoke test.

## Source version

- Public source commit: `ff5448abc6f9263d3a10f16a51e53e0d5c171fc4`
- Public commit title: `Initial public release of CCP-DP`
- Build identifier reported by the executable:
  `CCP-DP-v8.2-adaptive-routing-20260916`

The public source tree was clean before and after the reproduction.

## Full-dataset results

| Dataset | Length | Sequences | Exact | Success rate | Mean edit distance | Reconstruction rate | Algorithm time | Maximum RSS |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Srinivas-110 | 110 bp | 9,984 | 9,784 | 97.9968% | 0.049079 | 99.9554% | 2.03979 s | 198,708 KB |
| Chandak-108 | 108 bp | 1,466 | 1,398 | 95.3615% | 0.122101 | 99.8869% | 2.19428 s | 170,552 KB |
| Bar-Lev-140 | 140 bp | 10,000 | 9,874 | 98.7400% | 0.185200 | 99.8677% | 2.71242 s | 126,892 KB |
| Chandak-150 | 150 bp | 11,710 | 11,505 | 98.2494% | 0.052007 | 99.9653% | 1.88374 s | 122,664 KB |
| Erlich-152 | 152 bp | 72,000 | 71,984 | 99.9778% | 0.000222 | 99.9999% | 8.48303 s | 354,240 KB |

The three primary paper results were reproduced exactly:

| Dataset | Expected exact count | Reproduced exact count | Status |
|---|---:|---:|---|
| Srinivas-110 | 9,784 | 9,784 | Exact match |
| Chandak-108 | 1,398 | 1,398 | Exact match |
| Bar-Lev-140 | 9,874 | 9,874 | Exact match |

The reproduced Srinivas-110 and Chandak-108 prediction files were also
byte-identical to the preserved formal prediction files.

Runtime values can vary slightly with processor load, thread scheduling,
compiler version, and hardware. Reconstruction outputs and accuracy metrics
are the primary deterministic reproduction targets.

## Output validation

For every dataset:

- the number of reconstructed sequences matched the number of centers;
- every non-empty output sequence had the requested target length;
- the evaluator completed without an error;
- generated files remained excluded from Git.

## Prediction hashes

SHA-256 hashes of the reproduced prediction files:

| Dataset | SHA-256 |
|---|---|
| Srinivas-110 | `0df87f382cf90f2a83b4c70c4f5ca9d7a1252ea78ac4dbfa2ef0a61f757de6ac` |
| Chandak-108 | `0256c3780f557cbb339bc90f443bfb9b517b081c566e7920bea2485deb97a9ad` |
| Bar-Lev-140 | `880dce946f0df13f32f5407d528dd48fdb42c53a6433bb35db561e519ee2d991` |
| Chandak-150 | `4d1107d30db8385e901c37d405bf6828eadca7b475321b0cb3f01dd004c07ac2` |
| Erlich-152 | `fc0ada823ae929bae92ad1d7c4f23bf50df25701dac506f3072d45b3e4c7e0db` |

The generated summary CSV had SHA-256:

`e97e6328b781061e8a3b74ec6f854b2bea9e7c6ad96e1e957f86d58408111e4d`

## Validation environment

- Operating system: Ubuntu 20.04, Linux 5.15, x86-64
- Compiler: GCC/G++ 9.4.0
- Python: 3.11.15
- Processor: Intel Xeon Gold 6130H at 2.10 GHz
- CPU topology: 2 sockets, 16 cores per socket, 2 threads per core
- CCP-DP threads: 32

## Reproduction commands

Build the main executable from the repository root:

    mkdir -p build

    g++ \
      -std=c++17 \
      -O3 \
      -fopenmp \
      -Iinclude \
      src/main.cpp \
      -o build/ccpdp

A complete reconstruction follows this command structure:

    ./build/ccpdp \
      -i data/processed/DATASET/Clusters.txt \
      -l TARGET_LENGTH \
      -s "====" \
      -o predictions.txt \
      --diag diag.csv \
      --guide-output guides.txt \
      --mode auto \
      --jobs 32

Evaluate the output with:

    python scripts/evaluate_answers_extended.py \
      -o predictions.txt \
      -a data/processed/DATASET/Centers.txt

See `data/README.md` for dataset acquisition and checksum information, and
see the experiment-specific README files under `experiments/` for the
remaining reproduction workflows.
