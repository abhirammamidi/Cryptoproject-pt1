"""Task 1.1 - Shanks' Baby-Step Giant-Step."""
import math
from common import new_stats, Timer, pow_cost


def solve(g, h, p, order=None, **kwargs):
    n = (p - 1) if order is None else order
    stats = new_stats()
    x = None
    with Timer() as t:
        m = math.isqrt(n)
        if m * m < n:
            m += 1                                   # m = ceil(sqrt(n))

        # Baby steps: table[g^j] = j  for j = 0 .. m-1
        table = {}
        cur = 1
        for j in range(m):
            table.setdefault(cur, j)                 # keep the smallest j
            cur = cur * g % p
            stats["group_ops"] += 1
        stats["memory_items"] = len(table)           # peak table size

        # factor = g^(-m), using the Fermat inverse g^(p-2) since p is prime
        factor = pow(pow(g, p - 2, p), m, p)
        stats["group_ops"] += pow_cost(p - 2) + pow_cost(m)

        # Giant steps: gamma = h * factor^i ; stop at the first table hit
        gamma = h
        for i in range(m):
            j = table.get(gamma)
            if j is not None:
                x = (i * m + j) % n
                break
            gamma = gamma * factor % p
            stats["group_ops"] += 1
    stats["time_sec"] = t.elapsed

    if x is None or pow(g, x, p) != h:
        raise ValueError("BSGS failed: h is not in <g> or the order is wrong")
    return x, stats
