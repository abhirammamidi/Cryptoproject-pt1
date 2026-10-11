# Fix list (what was wrong, what was done, how it was checked)

| Reported problem | Cause | Fix |
|---|---|---|
| Index calculus F B=100, G B=180, G B=220 and the G comparison: "Relation collection limit reached" | a fixed trial cap far below the 10^6-10^7 candidates these B need; every candidate also cost a full `pow` | no default cap (optional `max_trials` / `max_seconds`); each candidate costs ONE multiplication and two values are tested per step (g^q = -1 gives p-v = g^(k+q)) |
| (found while fixing) dependent relations | with g = 5 in the factor base, relation(k+1) = relation(k) + log g | random stride walk v *= g^s |
| Linear sieve G: "31 rows for 748 unknowns" | B=220, M=700 gives too few smooth T(c1,c2); dense elimination | M grows automatically (old relations kept), then B if needed; sparse elimination mod q; the run table shows B and M actually used |
| Parameters differed between algorithms/data file | B_ic=180/260, B_ls=220, M_ls=700 in the data file vs 500/2000/700/900 in the runner | `params.py` is the single source (assignment values); data-file hints only with `--use-instance-params`; the used parameters are printed in every row |
| Brute force D "not run" | by design (extrapolation) | rate measured on D-sized numbers (`bruteforce.calibrate`); years printed; status EXTRAPOLATED |
| Pohlig-Hellman: no scaling study; inner BSGS not counted | - | stats split into outer / inner / total; `exp_ph_scaling.py` (cost follows sqrt of the largest prime, p fixed) |
| Rho: repeated iteration counts | one fixed walk function: all seeds fall onto the same giant cycle | keyed 3-way partition, new key per seed and restart; test asserts distinct counts |
| Runner crashed on a data file without `tasks` / `factors` / `INSTANCES` | hard-coded layout | `data_loader.py` (8 layouts, derived order and factorisation, clear errors) |
| Hard-coded task list | - | required tasks always run; every other cell of the A-H table is run if cheap, else marked infeasible with a reason |

## Code-quality pass
* `run_all.py`: cells skipped only because of the time/memory budget are now **NOT RUN** (with the estimated time),
  and **INFEASIBLE** is kept for > 1 year or more memory than the machine has. Estimates use the measured
  multiplication speed. Fixed "1 prime factors". Added `--cell-seconds`, `--ram-fraction`, `--infeasible-years`, `--out`.
* `pollard_rho.py`: the time limit uses the shared `Limits` class (removed a helper that read a private attribute).
* `data_loader.py`: the container check is two plain branches instead of one long and/or expression; errors chain with `from`.
* `common.py`: deterministic Miller-Rabin (13 prime bases, proven below 3.3*10^24), extra bases above that.
* Uniform docstrings on all six `solve()` functions (arguments, returned statistics, group_ops rule, exceptions).
* Lint (`ruff --select F,E7,E9,E741,B`): clean. (B905 `zip(strict=)` is not used because it needs Python 3.10.)
* New tests: run_all.py end to end on a foreign data file, unreadable data file, primality vs a sieve and vs
  famous strong pseudoprimes. 29 tests in total.

## Verified when this package was built (on the AI's machine; re-run on yours)
* `python run_all.py` on the supplied assignment_data.py: all required tasks PASS, 0 failures, both tables printed.
* The same on 16 different randomly generated instance sets spread over 8 data-file layouts (tuples, variables,
  per-field variables, list of dicts, strings, lowercase labels, hook function, missing instances), exit status 0.
* Edge cases: g not a generator, h outside the subgroup, composite p, empty or syntactically broken data file,
  tiny primes: clean N/A / error messages, no tracebacks.
* `part1/test_part1.py`: 29 tests pass. `part1/experiments/exp_previously_failed.py`: 9/9 pass.
* `part2/test_hashes.py` against the supplied SHA_VECTORS: all sections pass.
No timing in this package is a result to quote: run the programs yourselves and use your own measurements.
