# Part 1 — Verified experiment record

This record was generated from the uploaded A–H dataset and the revised
assignment-aware runner. Timings are from one run in this environment and
will vary by machine and Python version.

| Instance | Algorithm | Result | Time (s) | Group operations | Memory items | Notes |
|---|---|---:|---:|---:|---:|---|
| A | Brute Force | Verified | 0.0168 | 149,722 | 1 | Required |
| A | BSGS | Verified | 0.00020 | 1,053 | 833 | Required |
| B | Brute Force | Verified | 0.9976 | 9,336,324 | 1 | Required |
| B | BSGS | Verified | 0.00116 | 6,848 | 4,878 | Required |
| C | BSGS | Verified | 0.3793 | 1,266,031 | 859,318 | Required |
| C | Pollard Rho | Verified | 1.0881 | 3,639,029 | 2 | Required |
| C2 | Pollard Rho | Verified | 3.6250 | 12,512,420 | 2 | Required |
| D | Pohlig–Hellman | Verified | 0.00215 | 9,790 | 1,172 | Required |
| E | Pohlig–Hellman | Verified | 0.00028 | 2,910 | 4 | Required |
| F | Index Calculus | Verified | 0.4702 | 4,973,715 | 209 | Required |
| G | Index Calculus | Verified | 9.4946 | 70,192,335 | 1,086 | Required |
| G | Linear Sieve | Verified | 1.0233 | 522,934 | 2,328 | Required |

All returned values were independently checked by the runner with
`pow(g, x, p) == h`. The local test secret is not passed into any solver.

## Linear Sieve G rank interpretation

The run reported 1,017 pivots for 1,026 columns (rank deficiency 9).
This is not a full-column-rank system. However, the implementation separately
checks that every factor-base prime logarithm needed by the final individual
logarithm stage is determined. It then verifies each recovered prime log and
the final DLP result. The returned G logarithm passed the independent DLP
check. The report must state the measured rank accurately; do not call the
matrix full rank.

## Scope and limitations

- H is marked bonus-only in the supplied data and was not run because this
  package does not implement a dedicated H bonus solver.
- Extra cross-algorithm comparisons are not included in this required-task
  table. Use `python run_all.py --include-extra` only for supplementary tests.
- These timings and counters are implementation/environment measurements,
  not universal complexity claims.
