"""Task 1.0 - Brute-force discrete logarithm (the baseline)."""
from common import new_stats, Timer


def solve(g, h, p, order=None, **kwargs):
    n = (p - 1) if order is None else order
    stats = new_stats()
    stats["memory_items"] = 1            # only the running value 'cur'
    x_found = None
    with Timer() as t:
        cur = 1                          # cur == g^x at the top of the loop
        for x in range(n):
            if cur == h:
                x_found = x
                break
            cur = cur * g % p
            stats["group_ops"] += 1
    stats["time_sec"] = t.elapsed
    if x_found is None:
        raise ValueError("no solution: h is not a power of g")
    return x_found, stats
