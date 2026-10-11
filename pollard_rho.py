"""Task 1.2 - Pollard's rho for discrete logarithms (Floyd cycle finding)."""
import math
import random
from common import new_stats, Timer, pow_cost, solve_linear_congruence


def _step(xi, a, b, g, h, p, n):
    """One step of the 3-way partitioned walk. Returns the new (x, a, b).
    Partition rule (assumed): S0 = {x mod 3 == 0}, S1 = {1}, S2 = {2}."""
    s = xi % 3
    if s == 0:                       # S0: multiply by g  -> a + 1
        return xi * g % p, (a + 1) % n, b
    if s == 1:                       # S1: square         -> 2a, 2b
        return xi * xi % p, (2 * a) % n, (2 * b) % n
    return xi * h % p, a, (b + 1) % n   # S2: multiply by h -> b + 1


def solve(g, h, p, order=None, seed=None, max_restarts=50,
          max_candidates=10**6, **kwargs):
    n = (p - 1) if order is None else order
    rng = random.Random(seed)
    stats = new_stats()
    stats.update(restarts=0, iterations=0, candidates_checked=0)
    stats["memory_items"] = 2        # tortoise (x,a,b) and hare (x,a,b)
    x_found = None

    with Timer() as t:
        for attempt in range(max_restarts + 1):
            # Fresh random start x0 = g^a0 * h^b0
            a0, b0 = rng.randrange(n), rng.randrange(n)
            x0 = pow(g, a0, p) * pow(h, b0, p) % p
            stats["group_ops"] += pow_cost(a0) + pow_cost(b0) + 1

            tx, ta, tb = x0, a0, b0            # tortoise
            hx, ha, hb = x0, a0, b0            # hare
            limit = 4 * math.isqrt(n) + 100    # safety cap per attempt
            collided = False
            for _ in range(limit):
                tx, ta, tb = _step(tx, ta, tb, g, h, p, n)
                hx, ha, hb = _step(hx, ha, hb, g, h, p, n)
                hx, ha, hb = _step(hx, ha, hb, g, h, p, n)
                stats["group_ops"] += 3
                stats["iterations"] += 1
                if tx == hx:
                    collided = True
                    break
            if not collided:
                stats["restarts"] += 1
                continue

            # g^ta h^tb = g^ha h^hb  =>  x * (hb - tb) = (ta - ha)  (mod n)
            r = (hb - tb) % n
            s = (ta - ha) % n
            d = math.gcd(r, n)
            if r == 0 or d > max_candidates:   # d = n, or too many to test
                stats["restarts"] += 1
                continue
            for cand in solve_linear_congruence(r, s, n):
                stats["candidates_checked"] += 1
                stats["group_ops"] += pow_cost(cand)
                if pow(g, cand, p) == h:
                    x_found = cand
                    break
            if x_found is not None:
                break
            stats["restarts"] += 1             # no candidate verified
    stats["time_sec"] = t.elapsed

    if x_found is None:
        raise ValueError("Pollard rho failed after all restarts")
    return x_found, stats
