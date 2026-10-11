PART 1 — ENGINEERED WORKING PACKAGE
====================================

Requirements
------------
- Python 3.9+ (standard library only)
- No third-party discrete-log or cryptography packages are used.

Files
-----
assignment_data.py       Supplied A–H public parameters and local test secrets.
common.py                Shared arithmetic, factorisation, timing and stats helpers.
bruteforce.py             Brute-force DLP.
bsgs.py                   Baby-step giant-step.
pollard_rho.py            Pollard rho DLP.
pohlig_hellman.py         Pohlig–Hellman with prime-power lifting and CRT.
index_calculus.py         Index calculus for safe-prime groups.
linear_sieve.py           Linear-sieve index calculus.
run_all.py                Assignment-task-aware runner with independent verification.
test_part1.py             Fast smoke/unit tests.
experiments/results/      Generated CSV and Markdown run summaries.

Run
---
    python assignment_data.py
    python -m unittest -v test_part1.py
    python run_all.py

The default runner follows each instance's `tasks` field. It runs:
A: brute force, BSGS
B: brute force, BSGS
C: BSGS, Pollard rho
C2: Pollard rho
D/E: Pohlig–Hellman
F: index calculus
G: index calculus, linear sieve
H: bonus-only; a separate checker verifies the supplied local test secret, but no algorithmic H solver is claimed.

Optional:
    python run_all.py --only bsgs,rho
    python run_all.py --seed 123
    python run_all.py --include-extra

`--include-extra` is for exploratory comparisons only. It deliberately runs
selected solvers beyond the assignment's required task mapping.

Important algorithm notes
-------------------------
1. Every returned x in run_all.py is checked using:
       pow(g, x, p) == h
   A failed check is recorded as a failure and causes a non-zero exit status.

2. Linear Sieve G may report rank below the total number of columns. Some
   auxiliary H+c logarithm columns can remain free while all factor-base prime
   logarithms needed by Stage 3 are uniquely determined. The implementation
   explicitly checks every factor-base logarithm and verifies each recovered
   logarithm before using it. The final x is also checked by the runner.
   Do not describe this as a full-rank matrix in your report: quote the
   measured rank and explain which variables were actually required.

3. `group_ops` and `memory_items` are implementation-level estimates, not
   hardware-independent complexity measurements. State their counting
   conventions in your report. In particular, Python's `pow` calls and
   dictionary/list storage may not be fully represented by these counters.

4. Index-calculus and linear-sieve parameters are tuned for the supplied
   test instances. If your instructor provides different official parameters,
   use those instead and rerun the experiments.

5. H is a bonus instance. No H solver is claimed in this package.

Academic integrity
------------------
Use this as a code-review and testing aid. Read and understand each solver,
be prepared to explain it, and follow your course policy for acknowledging
friend/reference/AI assistance. Do not submit code you cannot explain.


Bonus H verification (known answer only)
----------------------------------------
    python bonus_h_check.py

This prints and verifies the H exponent stored in assignment_data.TEST_SECRETS.
It is verification data, not a discrete-log algorithm. Do not claim that
Pollard rho or another solver independently recovered H. A genuine solve of
the large prime-order component would require a substantial optimized
Pollard-rho computation and is not demonstrated by this package.


TERMINAL RESULTS TABLE
----------------------
Run:
    python run_all.py

The runner executes every algorithm assigned to each instance in
assignment_data.py and prints one compact table with status, recovered x,
runtime, group-operation count, and memory-item count. Successful answers
are independently checked with pow(g, x, p) == h. It also writes:
    experiments/results/results.csv
    experiments/results/summary.md

The default run intentionally follows the instance task lists. It does not
blindly run every method on every instance: some methods are incompatible
with certain group structures, while brute force or BSGS on large instances
can be impractical. Optional --include-extra forces extra experiments and
may take a long time or fail on unsuitable instances.

H is a bonus safe-prime instance. Its supplied test exponent is only checked
as a known-answer validation; it is not represented as an algorithmic DLP
solution. Do not claim H was solved by one of the algorithms.
