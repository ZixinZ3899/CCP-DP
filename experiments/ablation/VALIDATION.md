# Local validation

The ablation build compiled with `-O3 -march=native -std=c++17 -fopenmp
-Wall -Wextra -Wpedantic`.  A full one-repeat run on the supplied Srinivas-110
files produced the following accuracy metrics:

| Variant | Success (%) | Mean ED | Reconstruction (%) |
|---|---:|---:|---:|
| Guide only | 92.8586 | 0.180389 | 99.836010 |
| Full w/o Cross-LOO score | 94.1406 | 0.152444 | 99.861415 |
| Full w/o blind IDS likelihood | 97.9467 | 0.049880 | 99.954655 |
| Full with decoding of all clusters | 98.0068 | 0.048878 | 99.955565 |
| Local trellis only | 94.8618 | 0.136518 | 99.875892 |
| Full | 97.9968 | 0.049079 | 99.955383 |

The Full, Guide-only, w/o Cross-LOO, w/o blind-IDS, and decode-all values agree
with the previously reported v8.2/Srinivas behavior to rounding.  The local
trellis result is version-specific and must be reported from this rerun rather
than copied from the older table.

Runtime values from this local environment are intentionally not reported;
use the three-repeat server output in `ablation_summary.csv` for the paper.
