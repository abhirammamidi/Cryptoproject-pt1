"""Task 1.2 - Pollard's rho for discrete logarithms (Floyd cycle finding).

Walk: x_i = g^a_i * h^b_i with the assignment's three rules
    S0: x -> g*x   (a+1)        S1: x -> x^2  (2a, 2b)        S2: x -> h*x  (b+1)
WHICH rule applies is decided by a KEYED hash of x. Every seed (and every restart)
draws a fresh random key, so each run walks its own random mapping. With one fixed
mapping all seeds fall onto the same giant cycle and report identical iteration
counts, so repeated runs would not be independent samples.
"""
import math
import random
from common import new_stats, Timer, pow_cost, solve_linear_congruence, SolverFailure, Limits

_M64 = (1 << 64) - 1
_MULT = 0x9E3779B97F4A7C15            # odd 64-bit mixing constant (Fibonacci hashing)


def _make_step(g, h, p, n, salt):
    def step(x, a, b):
        y = ((x ^ (x >> 64) ^ salt) & _M64) * _MULT & _M64
        r = (y >> 32) % 3                                # keyed 3-way partition
        if r == 0:
            return x * g % p, (a + 1) % n, b             # S0
        if r == 1:
            return x * x % p, (2 * a) % n, (2 * b) % n   # S1
        return x * h % p, a, (b + 1) % n                 # S2
    return step


def solve(g, h, p, order=None, seed=None, max_restarts=100, max_candidates=10**6,
          max_seconds=None, **kwargs):
    """Pollard rho with Floyd cycle finding and the assignment's 3-way walk.

    Args:
        g, h, p: the problem g^x = h (mod p), p prime.
        order:   order n of g (default p - 1); the answer is returned modulo n.
        seed:    seed of this run (each seed and each restart uses its own random walk).
        max_restarts, max_candidates, max_seconds: failure bounds.
    Returns:
        (x, stats) with stats = {group_ops, time_sec, memory_items, iterations, restarts,
        candidates_checked}. group_ops = 3 per Floyd iteration + the starting exponentiations
        + the candidate checks; memory_items = 2 (tortoise and hare).
    Raises:
        SolverFailure if every restart fails or max_seconds is exhausted.
    """
    n = (p - 1) if order is None else order
    if n <= 1:
        raise ValueError("order must be > 1")
    rng = random.Random(seed)
    stats = new_stats()
    stats.update(restarts=0, iterations=0, candidates_checked=0)
    stats["memory_items"] = 2        # tortoise (x,a,b) and hare (x,a,b)
    x_found = None
    limit = 8 * math.isqrt(n) + 1000  # safety cap per attempt
    budget = Limits(max_seconds=max_seconds)

    with Timer() as t:
        for _attempt in range(max_restarts + 1):
            step = _make_step(g, h, p, n, rng.getrandbits(64))   # fresh random walk
            a0, b0 = rng.randrange(n), rng.randrange(n)
            x0 = pow(g, a0, p) * pow(h, b0, p) % p
            stats["group_ops"] += pow_cost(a0) + pow_cost(b0) + 1

            tx, ta, tb = hx, ha, hb = x0, a0, b0           # tortoise and hare start together
            collided = False
            for it in range(limit):
                tx, ta, tb = step(tx, ta, tb)               # tortoise: 1 step
                hx, ha, hb = step(*step(hx, ha, hb))        # hare: 2 steps
                if tx == hx:
                    collided = True
                    stats["iterations"] += it + 1
                    stats["group_ops"] += 3 * (it + 1)
                    break
                if (it & 0xFFFF) == 0 and budget.exceeded():
                    raise SolverFailure(f"rho stopped: {budget.exceeded()}", stats)
            if not collided:
                stats["iterations"] += limit
                stats["group_ops"] += 3 * limit
                stats["restarts"] += 1
                continue

            # g^ta h^tb = g^ha h^hb  =>  x * (hb - tb) = (ta - ha)  (mod n)
            r, s = (hb - tb) % n, (ta - ha) % n
            d = math.gcd(r, n)
            if r == 0 or d > max_candidates:   # d = n (no information) or too many candidates
                stats["restarts"] += 1
                continue
            for cand in solve_linear_congruence(r, s, n):   # d candidates (or none)
                stats["candidates_checked"] += 1
                stats["group_ops"] += pow_cost(cand)
                if pow(g, cand, p) == h:
                    x_found = cand
                    break
            if x_found is not None:
                break
            stats["restarts"] += 1             # no candidate verified: unusable collision
    stats["time_sec"] = t.elapsed

    if x_found is None:
        raise SolverFailure("Pollard rho failed after all restarts", stats)
    return x_found, stats

