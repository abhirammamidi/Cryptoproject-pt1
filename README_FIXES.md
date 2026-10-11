# Part 1 — reliability fixes

This package preserves the supplied assignment instances and solver interfaces, while adding bounded failure paths and more honest instrumentation.

## Fixes

- **Index calculus:** removed the no-op check; validates subgroup inputs; applies `max_trials` to initial relation collection, rank recovery, and individual-log relation collection. It now raises a diagnostic with relation/rank progress instead of potentially looping forever. This does not make an unsuitable smoothness bound magically produce smooth relations: increase `B` and/or the trial budget, and report failed parameter choices honestly.
- **Linear sieve:** adds a hard bound to individual-log relation search. The supplied G parameters (`B=220`, `M=700`) create 748 unknowns but only 31 usable structured rows in the reported run. That is a parameter/algorithm design shortfall, not a syntax error. The solver now fails clearly rather than hanging. A successful G run requires substantially more and independent relations plus a scalable sparse linear solver; simply increasing M with the current dense Gauss-Jordan routine can become prohibitively expensive.
- **Pohlig–Hellman:** now includes the internal BSGS operation count and peak baby-step table size in reported statistics.
- **Pollard rho:** uses a mixed residue for the walk partition to reduce dependence on the low bits of the group element. Continue using distinct seeds; report per-seed results and dispersion rather than treating identical iteration counts as proof of dependence.

## Run

From this directory:

```bash
python -m py_compile *.py experiments/*.py
python run_all.py
```

For index calculus, pass a realistic `max_trials` and sweep `B` independently. A failed parameter run is a valid experimental result; do not silently substitute another algorithm and label it index calculus.

## Important limitation

The current structured linear-sieve formulation still uses dense elimination. The original G configuration is underdetermined (31 rows for 748 columns). This package makes that limitation explicit and bounded; it does not claim that the original G linear-sieve instance is solved. A genuine successful implementation needs a larger structured-relation collection and sparse/modular elimination, or revised instructor-approved parameters.

The bundled experiment scripts are retained as parameter-sweep harnesses. Use `PYTHONPATH=..` when launching a script from the `experiments` directory. Their output should be included in the report with the exact B, M, seed, and trial budget used.
