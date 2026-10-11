"""Task 1.0 - Brute-force discrete logarithm (the baseline) + calibrated extrapolation."""
import time
from common import new_stats, Timer, SolverFailure


def solve(g, h, p, order=None, max_steps=None, **kwargs):
    """Brute force: try x = 0, 1, 2, ... until g^x = h.

    Args:
        g, h, p: the problem g^x = h (mod p), p prime.
        order:   order n of g (default p - 1); the answer is returned modulo n.
        max_steps: optional cap on the number of multiplications.
    Returns:
        (x, stats) with stats = {group_ops, time_sec, memory_items}.
        group_ops = x exactly (one multiplication per exponent tried); memory_items = 1.
    Raises:
        SolverFailure if max_steps is exhausted; ValueError if h is not a power of g.
    """
    n = (p - 1) if order is None else order
    stats = new_stats()
    stats["memory_items"] = 1            # only the running value 'cur'
    cap = n if max_steps is None else min(n, max_steps)
    x_found = None
    with Timer() as t:
        cur = 1                          # cur == g^x at the top of the loop
        for x in range(cap):
            if cur == h:
                x_found = x
                break
            cur = cur * g % p
            stats["group_ops"] += 1
    stats["time_sec"] = t.elapsed
    if x_found is None:
        if max_steps is not None and cap < n:
            raise SolverFailure(f"step budget {max_steps:,} exhausted before finding x", stats)
        raise ValueError("no solution: h is not a power of g")
    return x_found, stats


def calibrate(g, p, steps=300_000):
    """Measure seconds per multiplication g*cur mod p at the SIZE of this p
    (a bounded partial run: no solution is searched for)."""
    cur = 1
    t0 = time.perf_counter()
    for _ in range(steps):
        cur = cur * g % p
    return (time.perf_counter() - t0) / steps


SECONDS_PER_YEAR = 365.25 * 24 * 3600


def extrapolate(sec_per_op, n):
    """Return (average_seconds, worst_seconds, average_years, worst_years) for order n."""
    avg, worst = sec_per_op * n / 2, sec_per_op * n
    return avg, worst, avg / SECONDS_PER_YEAR, worst / SECONDS_PER_YEAR
