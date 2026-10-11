"""Task 1.5 - Linear sieve index calculus for a safe prime p = 2q + 1.

Relations come from (H + c1)(H + c2) = T(c1, c2)  (mod p)  with
H = isqrt(p) + 1, J = H^2 - p, T = J + (c1 + c2) H + c1 c2, 0 <= c1 <= c2 <= M.
Unknowns: log(small prime) for primes <= B, and log(H + c) for 0 <= c <= M.
"""
import math
import random
from common import new_stats, Timer, pow_cost, primes_up_to
from index_calculus import smooth_exponents, solve_mod_prime, lift_to_2q, individual_log


def _sieve_row(c1, H, J, M, primes, logp, power_cap):
    """For fixed c1, T(c1, c2) = A + c2 * D is linear in c2.
    Add log(l) at every c2 in [c1, M] where l^k divides T. Returns the array."""
    A = J + c1 * H
    D = H + c1
    width = M - c1 + 1
    acc = [0.0] * width
    for ell, lg in zip(primes, logp):
        if D % ell == 0:                  # T = A (mod l) for every c2
            if A % ell == 0:
                for i in range(width):
                    acc[i] += lg
            continue
        mod = ell
        while mod <= power_cap:
            r = (-A * pow(D, -1, mod)) % mod      # c2 = r (mod l^k)
            start = c1 + ((r - c1) % mod)
            for c2 in range(start, M + 1, mod):
                acc[c2 - c1] += lg
            mod *= ell
    return A, D, acc


def solve(g, h, p, order=None, B=700, M=900, seed=None, slack=1.0,
          power_cap=None, **kwargs):
    n = (p - 1) if order is None else order
    q = (p - 1) // 2
    if n != 2 * q:
        raise ValueError("this solver expects a safe prime and a generator g")
    primes = primes_up_to(B)
    if g not in primes:
        raise ValueError("anchor log_g(g) = 1 needs g in the factor base")
    L = len(primes)
    P = math.prod(primes)
    logp = [math.log(x) for x in primes]
    H = math.isqrt(p) + 1
    J = H * H - p
    power_cap = power_cap or (2 * H * (M + 1) + J)   # covers every T value
    ncols = L + M + 1                      # columns: primes, then H + c
    rng = random.Random(seed)
    stats = new_stats()
    stats.update(H=H, J=J, factor_base_size=ncols, candidates=0, sieve_hits=0)

    with Timer() as total:
        # ---------------- Stage 1: sieve for smooth T(c1, c2)
        rows, rhs = [], []
        with Timer() as t1:
            for c1 in range(M + 1):
                A, D, acc = _sieve_row(c1, H, J, M, primes, logp, power_cap)
                for i, val in enumerate(acc):
                    c2 = c1 + i
                    T = A + c2 * D
                    stats["candidates"] += 1
                    if val < math.log(T) - slack:
                        continue
                    stats["sieve_hits"] += 1
                    exps = smooth_exponents(T, primes, P)
                    if exps is None:
                        continue
                    row = {j: e for j, e in enumerate(exps) if e}
                    row[L + c1] = row.get(L + c1, 0) - 1     # - log(H + c1)
                    row[L + c2] = row.get(L + c2, 0) - 1     # - log(H + c2)
                    rows.append(row)
                    rhs.append(0)
            rows.append({primes.index(g): 1})                # anchor: log_g(g) = 1
            rhs.append(1)
        stats["relations"] = len(rows) - 1
        stats["memory_items"] = len(rows)

        # ---------------- Stage 2: linear algebra mod q, parity lift, verify
        with Timer() as t2:
            col_count = [0] * ncols
            for row in rows:
                for c in row:
                    col_count[c] += 1
            Xq, rank = solve_mod_prime(rows, rhs, ncols, q, col_count)
            stats["rank"] = rank
            stats["total_columns"] = ncols
            stats["rank_deficiency"] = ncols - rank
            stats["full_column_rank"] = (rank == ncols)
            # The auxiliary log_g(H+c) columns need not all be unique.
            # Stage 3 only requires the factor-base prime logarithms below.
            stats["factor_base_logs_determined"] = all(
                Xq[i] is not None for i in range(L)
            )
            logs = []
            for i, ell in enumerate(primes):
                if Xq[i] is None:
                    raise ValueError(f"log of {ell} undetermined: increase M or B")
                lg = lift_to_2q(Xq[i], ell, g, p, q)
                stats["group_ops"] += 2 * pow_cost(q)
                if pow(g, lg, p) != ell:
                    raise ValueError(f"wrong log for {ell}: check the relations")
                logs.append(lg)

        # ---------------- Stage 3: individual logarithm over the small primes
        with Timer() as t3:
            x = individual_log(h, g, p, n, primes, P, logs, rng, stats)
    stats.update(stage1_sec=t1.elapsed, stage2_sec=t2.elapsed,
                 stage3_sec=t3.elapsed, time_sec=total.elapsed)
    return x, stats
