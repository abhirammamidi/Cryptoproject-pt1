PART 1 - DLP REFERENCE PROJECT
===============================

This folder contains the Part 1 reference implementation and the assignment
instance data reorganised from the supplied group-member test-data file.

Files
-----
assignment_data.py   Canonical A-H data + separated TEST_SECRETS + Part-2 vectors
common.py             Shared number-theory helpers
bruteforce.py         Sequential DLP search
bsgs.py               Baby-step giant-step
pollard_rho.py        Pollard rho with restarts and collision handling
pohlig_hellman.py     Prime-power lifting + CRT
index_calculus.py     Index-calculus reference implementation
linear_sieve.py       Linear-sieve reference implementation
run_all.py            End-to-end runner

Why assignment_data.py is structured this way
----------------------------------------------
The numeric A-H values are retained from the supplied data, but the file is
engineered around the solver interface used by this project:

  INSTANCES[label]["p"]
  INSTANCES[label]["g"]
  INSTANCES[label]["h"]
  INSTANCES[label]["order"]
  INSTANCES[label]["factors"]

Algorithm-specific tuning is also attached where useful:
  B_ic, B_ls, M_ls

The secret exponents are NOT stored inside INSTANCES. They are isolated in
TEST_SECRETS so the solvers cannot accidentally consume the answer while
running. Running `python assignment_data.py` validates all DLP instances and
checks that the supplied test secrets reproduce h.

The Part-2 messages and SHA-256/SHA-512/SHA3-256 vectors are also retained.

IMPORTANT
---------
The uploaded assignment document says the official A-H numeric inputs are
supplied separately. This file is therefore a reference/test dataset based on
the supplied group-member data, not a claim that these are the official course
instances. Replace it with the instructor-provided data if different.

ACADEMIC INTEGRITY
------------------
The assignment requires submitted code to be understood and written by the
group, and asks students to disclose AI/reference implementations used. Treat
these files as a reference implementation: review, understand, test, and
rewrite/adapt where required before submission.
