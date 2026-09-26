# Third-party software notices

The root MIT License applies to the original CCP-DP implementation,
repository-authored adapters, evaluation tools, experiment orchestration, and
documentation, except where a file or directory states otherwise.

Third-party components retain their original copyright and license terms.
The root MIT License does not replace or override those terms.

## BBS

- Upstream: <https://github.com/GZHoffie/bbs>
- Pinned revision: `27c39ec5b631498b3c55073038fa6a59821f2974`
- License: MIT
- Local license: `baselines/bbs/LICENSE`

The repository contains the pinned upstream Rust source.

## CPL

- Upstream: <https://github.com/itaiorr/Deep-DNA-based-storage>
- Pinned revision: `167a68261719ea458e79c40fcc485ae768c1a64a`
- License: MIT
- Local license: `baselines/cpl/LICENSE`

The repository contains the required upstream C++ source and a local
benchmark adapter.

## TrellisBMA

- Pinned revision: `77cb3d3965548e0ff63a37ff5de13f96fa15f5a8`
- License: MIT
- Local license: `baselines/trellisbma/LICENSE`
- Attribution: `baselines/trellisbma/README_UPSTREAM.md`

The repository contains the upstream modules required by the local cluster
runner.

## MUSCLE

- Upstream: <https://github.com/rcedgar/muscle>
- Version used: MUSCLE 5.3
- License: GNU General Public License v3
- Local license copy: `baselines/muscle/LICENSE`

The MUSCLE executable and source code are not redistributed in this
repository. Users install MUSCLE separately using
`baselines/muscle/environment.yml`. The repository-authored wrapper invokes
the external executable as a separate process.

## ITR

- Upstream: <https://github.com/omersabary/Reconstruction>
- Pinned revision: `c50dec739bd2c7f18ac7678d921eecee02e86c6a`
- Upstream license statement: `License TBA`

The upstream ITR source is not redistributed. The script
`baselines/itr/fetch_upstream.sh` downloads the pinned revision into the
ignored `build/third_party/` directory. Users remain responsible for
complying with the upstream authors' terms.

## Deep-learning baselines

TReconLM, DNAFormer, and RobuSeqNet source trees and model weights are not
redistributed in this repository. The source-fetching instructions place them
under the ignored `build/` directory.

Repository-authored adapters, checksum records, split indices, and preserved
evaluation outputs remain available for reproducibility. Each upstream
project remains subject to its own license and usage terms.

## Datasets

Dataset files are not redistributed in this repository. Download and
preparation instructions are provided under `data/`. Dataset users must
follow the terms and citation requirements of the original data providers.

## Additional details

Exact revisions, checksums, and toolchain information are recorded in:

- `baselines/VERSIONS.md`
- `baselines/README.md`
- the README and license file inside each baseline directory
