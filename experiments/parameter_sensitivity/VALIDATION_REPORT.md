# Validation report

## Implementation checks

- The sensitivity executable compiles with `g++ -O3 -march=native -std=c++17 -fopenmp`.
- With no sensitivity flags, the Srinivas-110 and Chandak-108 result files are byte-identical to the production CCP-DP outputs.
- The runner's Myers Levenshtein implementation matched a conventional dynamic-programming implementation on 500 random sequence pairs.
- Python syntax checks and shell syntax checks passed.
- A two-dataset smoke test (100 clusters per dataset) completed successfully.
- The complete one-repeat, 36-run local validation completed and produced CSV, SVG, PDF, and PNG outputs.

## Default-result regression

| Dataset | Exact reconstruction | Mean ED | Structured decoding |
|---|---:|---:|---:|
| Srinivas-110 | 97.9968% | 0.0491 | 10.0% |
| Chandak-108 | 95.3615% | 0.1221 | 100.0% |

These are regression checks for the existing datasets, not new independent test results.

## Scope

The bundled preview used one repeat to verify the full analysis path. Use the server command in `README_CN.md` with three repeats for manuscript runtime statistics. The server configuration now contains the supplied absolute paths for Srinivas-110, Chandak-108, and Erlich-152. Erlich-152 was not available in the local workspace, so its default-result equivalence must be checked by the server smoke test before the three-dataset run.
