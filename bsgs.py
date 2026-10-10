"""Task 1.1 - Shanks' Baby-Step Giant-Step."""

import math
import time


def solve(g, h, p, order=None, **kwargs):
    n = p - 1 if order is None else order
    if n <= 0:
        raise ValueError("order must be positive")

    g %= p
    h %= p
    stats = {"group_ops": 0, "time_sec": 0.0, "memory_items": 0}

    start = time.perf_counter()

    if h == 1:
        stats["time_sec"] = time.perf_counter() - start
        return 0, stats

    m = math.isqrt(n)
    if m * m < n:
        m += 1

    # Baby table: g^j -> j.
    table = {}
    cur = 1
    for j in range(m):
        table.setdefault(cur, j)
        cur = cur * g % p
        stats["group_ops"] += 1

    stats["memory_items"] = len(table)

    # g^(-m) = (g^(-1))^m.
    g_inv = pow(g, -1, p)
    factor = pow(g_inv, m, p)
    stats["group_ops"] += max(0, g.bit_length() - 1) + max(0, m.bit_length() - 1)

    gamma = h
    for i in range(m + 1):
        j = table.get(gamma)
        if j is not None:
            candidate = i * m + j
            if candidate < n and pow(g, candidate, p) == h:
                stats["group_ops"] += max(0, candidate.bit_length() - 1)
                stats["time_sec"] = time.perf_counter() - start
                return candidate, stats

        if i != m:
            gamma = gamma * factor % p
            stats["group_ops"] += 1

    stats["time_sec"] = time.perf_counter() - start
    return None, stats


if __name__ == "__main__":
    x, st = solve(5, 8, 23, 22)
    print("x =", x)
    print("verification =", pow(5, x, 23) == 8)
    print(st)
