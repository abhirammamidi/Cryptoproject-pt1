"""Task 1.5 - Structured linear-sieve relation collection.

Relations use
(H+c1)(H+c2) = T(c1,c2) (mod p),
where H=floor(sqrt(p))+1 and J=H^2-p.
"""

import math
import random
import time
from common import primes_up_to
from index_calculus import _gauss_mod_prime, _lift_parity, _factor_over_base


def _sieve_row(c1, H, J, M, primes):
    A = J + c1*H
    D = H + c1
    width = M - c1 + 1
    scores = [0.0] * width

    for q in primes:
        if D % q == 0:
            if A % q == 0:
                for i in range(width):
                    scores[i] += math.log(q)
            continue

        for power in (q,):
            if power > A + (M+1)*D:
                break
            r = (-A * pow(D, -1, power)) % power
            first = c1 + ((r-c1) % power)
            for c2 in range(first, M+1, power):
                scores[c2-c1] += math.log(q)

    return A, D, scores


def solve(g, h, p, order=None, B=700, M=900, seed=None, slack=1.0, **kwargs):
    n = p - 1 if order is None else order
    q = (p - 1)//2
    if n != 2*q:
        raise ValueError("linear sieve implementation expects a safe prime and full order")
    base = primes_up_to(B)
    if g not in base:
        raise ValueError("g must be in the factor base")

    H = math.isqrt(p) + 1
    J = H*H - p
    rng = random.Random(seed)
    start = time.perf_counter()

    # Columns: small primes, then H+c.
    L = len(base)
    ncols = L + M + 1
    rows = []
    rhs = []

    for c1 in range(M+1):
        A, D, scores = _sieve_row(c1, H, J, M, base)
        for offset, score in enumerate(scores):
            c2 = c1 + offset
            T = A + c2*D
            if score + slack < math.log(T):
                continue

            exps = _factor_over_base(T, base)
            if exps is None:
                continue

            row = {i: e for i, e in enumerate(exps) if e}
            row[L+c1] = row.get(L+c1, 0) - 1
            row[L+c2] = row.get(L+c2, 0) - 1
            rows.append(row)
            rhs.append(0)

    # Anchor log_g(g)=1.
    rows.append({base.index(g): 1})
    rhs.append(1)

    if len(rows) < ncols:
        raise RuntimeError(
            f"insufficient independent structured relations: {len(rows)} rows for {ncols} unknowns; "
            "increase B or M"
        )

    dense_rows = []
    for row in rows:
        dense = [0]*ncols
        for c, v in row.items():
            dense[c] = v
        dense_rows.append(dense)

    logs_all, rank = _gauss_mod_prime(dense_rows, rhs, ncols, q)
    if logs_all is None:
        raise RuntimeError("linear-sieve relation matrix is rank deficient; increase B/M")

    logs = []
    for i, prime in enumerate(base):
        logs.append(_lift_parity(logs_all[i], prime, g, p, q))

    # Individual logarithm using the recovered factor-base logs.
    tries = 0
    while True:
        k = rng.randrange(n)
        value = h * pow(g, k, p) % p
        tries += 1
        exps = _factor_over_base(value, base)
        if exps is None:
            continue
        x = (sum(e*l for e, l in zip(exps, logs)) - k) % n
        if pow(g, x, p) == h:
            break

    total = time.perf_counter() - start
    stats = {
        "group_ops": tries,
        "time_sec": total,
        "memory_items": len(rows),
        "factor_base_size": ncols,
        "relations": len(rows)-1,
        "rank": rank,
        "stage1_sec": total,
        "stage2_sec": 0.0,
        "stage3_sec": 0.0,
        "stage3_tries": tries,
        "H": H,
        "J": J,
    }
    return x, stats
