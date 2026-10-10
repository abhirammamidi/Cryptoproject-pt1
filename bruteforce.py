"""Task 1.0 - trial multiplication discrete logarithm."""

import time


def solve(g, h, p, order=None, **kwargs):
    if p <= 1:
        raise ValueError("p must be > 1")
    n = p - 1 if order is None else order
    if n <= 0:
        raise ValueError("order must be positive")

    g %= p
    h %= p
    cur = 1
    stats = {"group_ops": 0, "time_sec": 0.0, "memory_items": 1}

    start = time.perf_counter()
    for x in range(n):
        if cur == h:
            stats["time_sec"] = time.perf_counter() - start
            return x, stats
        cur = cur * g % p
        stats["group_ops"] += 1

    stats["time_sec"] = time.perf_counter() - start
    return None, stats


if __name__ == "__main__":
    x, st = solve(5, 8, 23, 22)
    print("x =", x)
    print("verification =", pow(5, x, 23) == 8)
    print(st)
