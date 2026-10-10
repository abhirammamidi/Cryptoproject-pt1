"""Task 1.4 - Index calculus for safe-prime groups.

This implementation follows the assignment's safe-prime guidance:
solve the relation system modulo q=(p-1)/2, then recover the missing
parity by testing x and x+q.
"""

import math
import random
import time
from common import primes_up_to


def _factor_over_base(value, base):
    exponents = [0] * len(base)
    x = value
    for i, q in enumerate(base):
        while x % q == 0:
            x //= q
            exponents[i] += 1
        if x == 1:
            return exponents
    return exponents if x == 1 else None


def _gauss_mod_prime(rows, rhs, cols, q):
    """Dense Gauss-Jordan elimination over F_q."""
    a = [[v % q for v in row] + [b % q] for row, b in zip(rows, rhs)]
    rank = 0
    pivots = []

    for col in range(cols):
        pivot = None
        for r in range(rank, len(a)):
            if a[r][col]:
                pivot = r
                break
        if pivot is None:
            continue

        a[rank], a[pivot] = a[pivot], a[rank]
        inv = pow(a[rank][col], q - 2, q)
        a[rank] = [(v * inv) % q for v in a[rank]]

        for r in range(len(a)):
            if r == rank or a[r][col] == 0:
                continue
            f = a[r][col]
            a[r] = [(u - f*v) % q for u, v in zip(a[r], a[rank])]

        pivots.append(col)
        rank += 1
        if rank == len(a):
            break

    if rank < cols:
        return None, rank

    solution = [0] * cols
    for r, col in enumerate(pivots[:cols]):
        solution[col] = a[r][-1]
    return solution, rank


def _lift_parity(xq, element, g, p, q):
    for candidate in (xq, xq + q):
        if pow(g, candidate, p) == element:
            return candidate
    raise ValueError("parity lift failed")


def solve(g, h, p, order=None, B=500, relation_factor=2.0,
          seed=None, max_trials=None, **kwargs):
    n = p - 1 if order is None else order
    q = (p - 1) // 2

    if n != 2*q:
        raise ValueError("index_calculus implementation expects order p-1 for a safe prime")
    if not (pow(2, 1, p) and all(True for _ in [0])):
        pass

    base = primes_up_to(B)
    L = len(base)
    if g not in base:
        raise ValueError("g must be in the factor base for this implementation")

    rng = random.Random(seed)
    start = time.perf_counter()
    rows = []
    rhs = []
    trials = 0
    target = max(L + 5, int(math.ceil(relation_factor * L)))

    while len(rows) < target:
        if max_trials is not None and trials >= max_trials:
            raise RuntimeError("relation collection limit reached")
        k = rng.randrange(1, n)
        value = pow(g, k, p)
        trials += 1
        exps = _factor_over_base(value, base)
        if exps is not None:
            rows.append(exps)
            rhs.append(k % q)

    stage1 = time.perf_counter() - start

    t2 = time.perf_counter()
    logs_q, rank = _gauss_mod_prime(rows, rhs, L, q)
    if logs_q is None:
        # Continue collecting until the system reaches full rank.
        while rank < L:
            k = rng.randrange(1, n)
            value = pow(g, k, p)
            trials += 1
            exps = _factor_over_base(value, base)
            if exps is not None:
                rows.append(exps)
                rhs.append(k % q)
            logs_q, rank = _gauss_mod_prime(rows, rhs, L, q)
    logs = [_lift_parity(v, prime, g, p, q) for v, prime in zip(logs_q, base)]
    stage2 = time.perf_counter() - t2

    # Individual logarithm: find k for which h*g^k is B-smooth.
    t3 = time.perf_counter()
    tries = 0
    while True:
        k = rng.randrange(0, n)
        value = h * pow(g, k, p) % p
        tries += 1
        exps = _factor_over_base(value, base)
        if exps is None:
            continue

        x = (sum(e*l for e, l in zip(exps, logs)) - k) % n
        if pow(g, x, p) == h:
            break

    stage3 = time.perf_counter() - t3
    total = time.perf_counter() - start

    stats = {
        "group_ops": trials + tries,
        "time_sec": total,
        "memory_items": len(rows),
        "factor_base_size": L,
        "relations": len(rows),
        "rank": rank,
        "stage1_sec": stage1,
        "stage2_sec": stage2,
        "stage3_sec": stage3,
        "stage3_tries": tries,
    }
    return x, stats
