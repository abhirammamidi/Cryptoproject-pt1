"""Task 1.4 - Index calculus for p = 2q + 1 (safe prime), g a generator.

Speed-ups over the textbook loop (all exact, none changes the algorithm):
  * strided exponents: v = g^(k+s) = v*u (u = g^s, s random) costs ONE multiplication,
    not a full pow. (A stride of 1 would be wrong when g itself is in the factor base:
    relation(k+1) = relation(k) + log g, i.e. linearly dependent rows.)
  * the sign trick: g^q = -1, so p - v = g^(k+q): every step tests TWO values for free.
Failures are controlled: a budget or rank problem raises IndexCalculusFailure with
diagnostics instead of looping forever.
"""
import math
import random
import time
from common import new_stats, Timer, pow_cost, primes_up_to, SolverFailure, Limits

_WINDOW = 2048


class IndexCalculusFailure(SolverFailure):
    pass


# ------------------------------------------------------------ smoothness
def smooth_exponents(v, primes, P):
    """Exponent list of v over `primes` if v is smooth, else None.
    P = product of all primes; gcd(v, P) rejects non-smooth v cheaply."""
    d = math.gcd(v, P)
    if d == 1:
        return None if v != 1 else [0] * len(primes)
    w = v
    while d > 1:                         # strip every base prime dividing v
        w //= d
        d = math.gcd(w, d)
    if w != 1:
        return None                      # a prime > B remains
    exps = []
    for q in primes:                     # v is smooth: trial-divide for exponents
        e = 0
        while v % q == 0:
            v //= q
            e += 1
        exps.append(e)
    return exps


# ------------------------------------------------------------ linear algebra
def solve_mod_prime(rows, rhs, ncols, q, col_count=None):
    """Solve  sum(row[c] * X[c]) = rhs  (mod prime q)  by sparse Gauss-Jordan.
    rows are dicts {column: coefficient}. Returns (X, rank), where X[c] is
    None when column c is not uniquely determined (a rank deficiency).
    Invariant: every stored pivot row has coefficient 1 on its own pivot and
    0 on every other pivot column, so one pass reduces any new row."""
    pivots = {}                                  # pivot col -> [row dict, rhs]
    for row, b in zip(rows, rhs):
        r = {c: v % q for c, v in row.items() if v % q}
        b %= q
        for c in [c for c in r if c in pivots]:  # reduce by existing pivots
            f = r.get(c)
            if not f:
                continue
            prow, pb = pivots[c]
            for cc, vv in prow.items():
                nv = (r.get(cc, 0) - f * vv) % q
                if nv:
                    r[cc] = nv
                else:
                    r.pop(cc, None)
            b = (b - f * pb) % q
        if not r:
            continue                             # row was a combination of others
        if col_count is None:
            c = min(r)
        else:                                    # pivot on the rarest column
            c = min(r, key=lambda col: (col_count[col], col))
        inv = pow(r[c], q - 2, q)                # Fermat inverse, q prime
        r = {cc: vv * inv % q for cc, vv in r.items()}
        b = b * inv % q
        for other in pivots.values():            # clear column c elsewhere
            f = other[0].get(c)
            if f:
                orow = other[0]
                for cc, vv in r.items():
                    nv = (orow.get(cc, 0) - f * vv) % q
                    if nv:
                        orow[cc] = nv
                    else:
                        orow.pop(cc, None)
                other[1] = (other[1] - f * b) % q
        pivots[c] = [r, b]

    X = [None] * ncols
    for c, (prow, pb) in pivots.items():
        if len(prow) == 1:                       # no free columns left in the row
            X[c] = pb
    return X, len(pivots)


# ------------------------------------------------------------ parity lift
def lift_to_2q(xq, element, g, p, q):
    """Given log_g(element) mod q, return it mod 2q using Euler's criterion:
    element^q == 1 (mod p) exactly when its log is even (g a generator)."""
    want_even = pow(element, q, p) == 1
    for cand in (xq, xq + q):
        if (cand % 2 == 0) == want_even:
            return cand
    raise AssertionError("unreachable")


# ------------------------------------------------------------ stage 1
def collect_relations(g, p, n, q, primes, P, need, rows, rhs, used, rng, stats, limits):
    """Append smooth relations  k = sum e_i log q_i  until len(rows) >= need."""
    while len(rows) < need:
        k0, s = rng.randrange(n), rng.randrange(2, n)
        v, u = pow(g, k0, p), pow(g, s, p)               # walk v = g^(k0 + s*j), j = 0, 1, ...
        stats["group_ops"] += pow_cost(k0) + pow_cost(s) + _WINDOW
        for j in range(_WINDOW):
            kj = k0 + s * j
            e = smooth_exponents(v, primes, P)           # value g^k
            if e is not None:
                k = kj % n
                if k not in used and any(e):
                    used.add(k)
                    rows.append({i: x for i, x in enumerate(e) if x})
                    rhs.append(k)
            e = smooth_exponents(p - v, primes, P)       # value -g^k = g^(k+q)
            if e is not None:
                k = (kj + q) % n
                if k not in used and any(e):
                    used.add(k)
                    rows.append({i: x for i, x in enumerate(e) if x})
                    rhs.append(k)
            v = v * u % p
        stats["candidates"] += 2 * _WINDOW
        why = limits.exceeded(stats["candidates"])
        if why:
            stats["relations"] = len(rows)
            raise IndexCalculusFailure(
                f"relation collection stopped: {why} (have {len(rows)}/{need} relations "
                f"after {stats['candidates']:,} candidates); raise B or the budget", stats)


# ------------------------------------------------------------ stage 3
def individual_log(h, g, p, n, q, primes, P, logs, rng, stats, limits):
    """Find k with h*g^k smooth; then x = sum(e_i*log q_i) - k (mod n)."""
    tested = 0
    while True:
        k0, s = rng.randrange(n), rng.randrange(2, n)
        v, u = h * pow(g, k0, p) % p, pow(g, s, p)
        stats["group_ops"] += pow_cost(k0) + pow_cost(s) + 1 + _WINDOW
        for j in range(_WINDOW):
            kj = k0 + s * j
            for val, k in ((v, kj), (p - v, kj + q)):
                e = smooth_exponents(val, primes, P)
                if e is not None:
                    x = (sum(ei * li for ei, li in zip(e, logs) if ei) - k) % n
                    if pow(g, x, p) == h:
                        stats["stage3_tries"] = tested + 2 * j + 2
                        return x
            v = v * u % p
        tested += 2 * _WINDOW
        why = limits.exceeded(stats["candidates"] + tested)
        if why:
            raise IndexCalculusFailure(f"individual-log search stopped: {why}", stats)


def suggest_B(p):
    """Heuristic starting bound B ~ exp(sqrt(0.35 ln p ln ln p)) (used only if no B is given)."""
    lp = math.log(p)
    return max(50, int(math.exp(math.sqrt(0.35 * lp * math.log(lp)))))


# ------------------------------------------------------------ the solver
def solve(g, h, p, order=None, B=500, relation_factor=2.0, seed=None,
          max_trials=None, max_seconds=None, **kwargs):
    """Index calculus for a safe prime p = 2q + 1 with a generator g (three stages).

    Args:
        g, h, p: the problem g^x = h (mod p), p prime.
        order:   order n of g (default p - 1); the answer is returned modulo n.
        B:       smoothness bound (factor base = primes <= B).
        relation_factor: relations collected per unknown before the first solve (default 2).
        seed, max_trials, max_seconds: reproducibility and failure bounds.
    Returns:
        (x, stats) with stats = {group_ops, time_sec, memory_items, factor_base_size, B,
        candidates, relations, rank, full_column_rank, factor_base_logs_determined,
        stage1_sec, stage2_sec, stage3_sec, stage3_tries}. group_ops = every modular
        multiplication of the walks + square-and-multiply cost of each pow;
        memory_items = number of stored relations.
    Raises:
        ValueError if p is not a safe prime or g is not a generator;
        IndexCalculusFailure if max_trials / max_seconds is exhausted.
    """
    n = (p - 1) if order is None else order
    q = (p - 1) // 2
    if n != 2 * q or 2 * q + 1 != p:
        raise ValueError("this solver expects a safe prime p = 2q + 1 and order p - 1")
    if pow(g, q, p) != p - 1:
        raise ValueError("g must be a generator (g^q must equal -1 mod p)")
    if not 1 <= h < p:
        raise ValueError("h out of range")
    rng = random.Random(seed)
    primes = primes_up_to(B)
    L = len(primes)
    if L < 2:
        raise ValueError("B too small")
    P = math.prod(primes)
    stats = new_stats()
    stats.update(factor_base_size=L, candidates=0, relations=0, B=B, rank=0)
    limits = Limits(max_trials, max_seconds)

    with Timer() as total:
        rows, rhs, used = [], [], set()
        target = int(relation_factor * L)
        stage1 = stage2 = 0.0
        while True:
            # ---------------- Stage 1: relation collection
            t_start = time.perf_counter()
            try:
                collect_relations(g, p, n, q, primes, P, target, rows, rhs, used,
                                  rng, stats, limits)
            finally:
                stage1 += time.perf_counter() - t_start
            # ---------------- Stage 2: linear algebra mod q, then parity lift
            with Timer() as t2:
                Xq, rank = solve_mod_prime(rows, rhs, L, q)
                ok = all(v is not None for v in Xq)
                if ok:
                    logs = [lift_to_2q(Xq[i], primes[i], g, p, q) for i in range(L)]
                    stats["group_ops"] += L * pow_cost(q)
            stage2 += t2.elapsed
            stats["rank"] = rank
            stats["relations"] = len(rows)
            if ok:
                break
            target = len(rows) + max(10, L // 5)     # rank deficient: collect more
        stats["memory_items"] = len(rows)
        stats["full_column_rank"] = rank >= L
        stats["factor_base_logs_determined"] = True
        # ---------------- Stage 3: individual logarithm
        with Timer() as t3:
            x = individual_log(h, g, p, n, q, primes, P, logs, rng, stats, limits)
    stats.update(stage1_sec=stage1, stage2_sec=stage2, stage3_sec=t3.elapsed,
                 time_sec=total.elapsed)
    return x, stats
