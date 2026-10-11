"""Task 1.1 - Shanks' Baby-Step Giant-Step."""
import math
from common import new_stats, Timer, pow_cost, SolverFailure


def solve(g, h, p, order=None, max_table=None, **kwargs):
    """Shanks' baby-step giant-step: x = i*m + j with m = ceil(sqrt(n)).

    Args:
        g, h, p: the problem g^x = h (mod p), p prime.
        order:   order n of g (default p - 1); the answer is returned modulo n.
        max_table: optional cap on the number of baby-step table entries.
    Returns:
        (x, stats) with stats = {group_ops, time_sec, memory_items}.
        group_ops = m baby-step multiplications + the cost of g^(-m) (square-and-multiply)
        + one multiplication per giant step; memory_items = peak table size (= m).
    Raises:
        SolverFailure if the table would exceed max_table; ValueError if no solution exists.
    """
    n = (p - 1) if order is None else order
    m = math.isqrt(n)
    if m * m < n:
        m += 1                                       # m = ceil(sqrt(n))
    if max_table is not None and m > max_table:
        raise SolverFailure(f"BSGS needs a table of {m:,} entries (limit {max_table:,})",
                            {"needed_entries": m})
    stats = new_stats()
    x = None
    with Timer() as t:
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
