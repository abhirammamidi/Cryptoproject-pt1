"""Task 1.4 - Index calculus for p = 2q + 1 (safe prime), g a generator."""
import math
import random
from common import new_stats, Timer, pow_cost, primes_up_to


# ------------------------------------------------------------ smoothness
def smooth_exponents(v, primes, P):
    """Exponent list of v over `primes` if v is smooth, else None.
    P = product of all primes; the gcd test rejects non-smooth v cheaply."""
    d = math.gcd(v, P)
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


def individual_log(h, g, p, n, primes, P, logs, rng, stats):
    """Stage 3: find k with h*g^k smooth; x = sum(e_i*log q_i) - k (mod n)."""
    tries = 0
    while True:
        k = rng.randrange(n)
        v = h * pow(g, k, p) % p
        stats["group_ops"] += pow_cost(k) + 1
        tries += 1
        exps = smooth_exponents(v, primes, P)
        if exps is None:
            continue
        x = (sum(e * L for e, L in zip(exps, logs) if e) - k) % n
        stats["group_ops"] += pow_cost(x)
        if pow(g, x, p) == h:
            stats["stage3_tries"] = tries
            return x


# ------------------------------------------------------------ the solver
def solve(g, h, p, order=None, B=500, relation_factor=2.0, seed=None, **kwargs):
    n = (p - 1) if order is None else order
    q = (p - 1) // 2
    if n != 2 * q:
        raise ValueError("this solver expects a safe prime and a generator g")
    rng = random.Random(seed)
    primes = primes_up_to(B)
    L = len(primes)
    P = math.prod(primes)
    stats = new_stats()
    stats.update(factor_base_size=L, candidates=0, relations=0)

    with Timer() as total:
        rows, rhs = [], []
        target = int(relation_factor * L)
        stage1 = stage2 = 0.0
        while True:
            # ---------------- Stage 1: relation collection
            with Timer() as t1:
                while len(rows) < target:
                    k = rng.randrange(1, n)
                    v = pow(g, k, p)
                    stats["group_ops"] += pow_cost(k)
                    stats["candidates"] += 1
                    exps = smooth_exponents(v, primes, P)
                    if exps is not None:
                        rows.append({i: e for i, e in enumerate(exps) if e})
                        rhs.append(k)
            stage1 += t1.elapsed
            # ---------------- Stage 2: linear algebra mod q, then parity lift
            with Timer() as t2:
                Xq, rank = solve_mod_prime(rows, rhs, L, q)
                ok = all(v is not None for v in Xq)
                if ok:
                    logs = [lift_to_2q(Xq[i], primes[i], g, p, q) for i in range(L)]
                    for i in range(L):
                        stats["group_ops"] += pow_cost(q)
            stage2 += t2.elapsed
            stats["rank"] = rank
            if ok:
                break
            target = len(rows) + max(10, L // 5)     # rank deficient: collect more
        stats["relations"] = len(rows)
        stats["memory_items"] = len(rows)
        # ---------------- Stage 3: individual logarithm
        with Timer() as t3:
            x = individual_log(h, g, p, n, primes, P, logs, rng, stats)
    stats.update(stage1_sec=stage1, stage2_sec=stage2, stage3_sec=t3.elapsed,
                 time_sec=total.elapsed)
    return x, stats
