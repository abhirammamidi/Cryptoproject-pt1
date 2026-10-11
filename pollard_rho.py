"""Task 1.2 - Pollard rho for discrete logarithms, Floyd cycle detection."""

import math
import random
import time


def _step(state, g, h, p, n):
    x, a, b = state
    # A mixed 3-way partition avoids relying on the low residue bits alone.
    z = (x ^ (x >> 17) ^ (x >> 31)) % 3
    r = z
    if r == 0:
        return (x * g) % p, (a + 1) % n, b
    if r == 1:
        return (x * x) % p, (2 * a) % n, (2 * b) % n
    return (x * h) % p, a, (b + 1) % n


def solve(g, h, p, order=None, seed=None, max_restarts=100, **kwargs):
    n = p - 1 if order is None else order
    if n <= 1:
        raise ValueError("order must be > 1")

    rng = random.Random(seed)
    stats = {
        "group_ops": 0, "time_sec": 0.0, "memory_items": 2,
        "restarts": 0, "iterations": 0, "candidates_checked": 0
    }

    start = time.perf_counter()

    for restart in range(max_restarts):
        a0 = rng.randrange(n)
        b0 = rng.randrange(n)
        x0 = pow(g, a0, p) * pow(h, b0, p) % p
        stats["group_ops"] += max(0, a0.bit_length()-1) + max(0, b0.bit_length()-1) + 1

        tortoise = (x0, a0, b0)
        hare = tortoise
        limit = max(1000, 8 * math.isqrt(n) + 100)

        for _ in range(limit):
            tortoise = _step(tortoise, g, h, p, n)
            hare = _step(_step(hare, g, h, p, n), g, h, p, n)
            stats["group_ops"] += 3
            stats["iterations"] += 1

            if tortoise[0] != hare[0]:
                continue

            _, a1, b1 = tortoise
            _, a2, b2 = hare

            numerator = (a1 - a2) % n
            denominator = (b2 - b1) % n
            d = math.gcd(denominator, n)

            if d == n:
                break

            # Solve denominator*x = numerator (mod n).
            if numerator % d:
                break

            n1 = n // d
            if n1 == 1:
                candidates = range(n)
            else:
                den1 = denominator // d
                num1 = numerator // d
                inv = pow(den1, -1, n1)
                base = (num1 * inv) % n1
                candidates = (base + k*n1 for k in range(d))

            for candidate in candidates:
                stats["candidates_checked"] += 1
                if pow(g, candidate, p) == h:
                    stats["group_ops"] += max(0, candidate.bit_length()-1)
                    stats["restarts"] = restart
                    stats["time_sec"] = time.perf_counter() - start
                    return candidate, stats

            break

        stats["restarts"] += 1

    stats["time_sec"] = time.perf_counter() - start
    return None, stats


if __name__ == "__main__":
    x, st = solve(5, 8, 23, 22, seed=7)
    print("x =", x)
    print("verification =", x is not None and pow(5, x, 23) == 8)
    print(st)
