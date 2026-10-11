"""Task 1.5 - Linear sieve index calculus for a safe prime p = 2q + 1.

Relations come from (H + c1)(H + c2) = T(c1, c2)  (mod p)  with
H = isqrt(p) + 1, J = H^2 - p, T = J + (c1 + c2) H + c1 c2, 0 <= c1 <= c2 <= M.
Unknowns: log(small prime) for primes <= B (columns 0..L-1) and log(H + c) for
0 <= c <= M (column L + c).  A homogeneous system, so the anchor log_g(g) = 1 is added.

If the requested (B, M) yields too few independent relations, the sieve range M is
GROWN (only the new pairs are sieved; old relations are kept) up to max_expand times,
and the final M is reported in the stats. Nothing is hidden: stats['M_requested'],
stats['M_final'] and stats['expansions'] say exactly what happened.
"""
import math
import random
import time
from common import new_stats, Timer, pow_cost, primes_up_to, SolverFailure, Limits
from index_calculus import smooth_exponents, solve_mod_prime, lift_to_2q, individual_log


class LinearSieveFailure(SolverFailure):
    pass


def _sieve_segment(c1, lo, hi, H, J, primes, logp):
    """For fixed c1, T(c1, c2) = A + D*c2 is LINEAR in c2. Return (A, D, acc) where acc[i]
    sums log(l) over every prime power l^k dividing T(c1, lo + i), for c2 in [lo, hi]."""
    A, D = J + c1 * H, H + c1
    width = hi - lo + 1
    acc = [0.0] * width
    Tmax = A + D * hi
    for ell, lg in zip(primes, logp):
        if D % ell == 0:                         # then T = A (mod l) for every c2
            if A % ell == 0:
                for i in range(width):
                    acc[i] += lg
            continue
        K, mod = 1, ell
        while mod * ell <= Tmax:                 # highest power that can divide some T
            mod *= ell
            K += 1
        t = (-A * pow(D, -1, mod)) % mod         # one inverse per prime, reused for every power
        m = ell
        for _ in range(K):
            start = lo + ((t % m - lo) % m)      # first c2 >= lo with c2 = -A/D (mod m)
            for c2 in range(start, hi + 1, m):
                acc[c2 - lo] += lg
            m *= ell
    return A, D, acc


def solve(g, h, p, order=None, B=700, M=900, seed=None, slack=1.0, auto_expand=True,
          max_expand=3, growth=1.35, max_B_bumps=8, B_growth=1.6, max_trials=None,
          max_seconds=None, **kwargs):
    """Linear-sieve index calculus for a safe prime p = 2q + 1 with a generator g.

    Args:
        g, h, p: the problem g^x = h (mod p), p prime.
        order:   order n of g (default p - 1); the answer is returned modulo n.
        B, M:    factor-base bound and sieve range 0 <= c1 <= c2 <= M.
        auto_expand: if (B, M) give too few relations, grow M (max_expand times) and then,
                 as a last resort, raise B and restart (max_B_bumps times).
        seed, max_trials, max_seconds: reproducibility and failure bounds.
    Returns:
        (x, stats) with stats = {group_ops, time_sec, memory_items, H, J, B_requested, B_final,
        M_requested, M_final, expansions, B_bumps, factor_base_size (unknowns), candidates,
        sieve_hits, relations, rank, full_column_rank, factor_base_logs_determined,
        stage1_sec, stage2_sec, stage3_sec}. memory_items = number of stored relations.
    Raises:
        ValueError if p is not a safe prime, g is not a generator, or g is not in the factor base;
        LinearSieveFailure if no parameter growth helps or max_seconds is exhausted.
    """
    limits = Limits(max_trials, max_seconds)
    B0, M0, last = B, M, None
    for bump in range(max_B_bumps + 1):
        try:
            x, st = _solve_once(g, h, p, order, B, M, seed, slack, auto_expand, max_expand, growth, limits)
            st["B_requested"], st["B_final"], st["B_bumps"] = B0, B, bump
            st["M_requested"] = M0                          # what the caller asked for (M_final = what was used)
            return x, st
        except LinearSieveFailure as exc:
            last = exc
            if not auto_expand or bump == max_B_bumps or not getattr(exc, "retryable", False):
                raise
            B = max(B + 1, int(B * B_growth))
            M = max(M, int(exc.stats.get("M_final", M)))      # keep the grown sieve range as well
    raise last


def _solve_once(g, h, p, order, B, M, seed, slack, auto_expand, max_expand, growth, limits):
    n = (p - 1) if order is None else order
    q = (p - 1) // 2
    if n != 2 * q or 2 * q + 1 != p:
        raise ValueError("this solver expects a safe prime p = 2q + 1 and order p - 1")
    if pow(g, q, p) != p - 1:
        raise ValueError("g must be a generator (g^q must equal -1 mod p)")
    primes = primes_up_to(B)
    index = {ell: i for i, ell in enumerate(primes)}
    if g not in index:
        raise ValueError("anchor log_g(g) = 1 needs g in the factor base (g prime and <= B)")
    L = len(primes)
    P = math.prod(primes)
    logp = [math.log(x) for x in primes]
    H = math.isqrt(p) + 1
    J = H * H - p
    rng = random.Random(seed)
    stats = new_stats()
    stats.update(H=H, J=J, B=B, M_requested=M, M_final=M, expansions=0, candidates=0,
                 sieve_hits=0, relations=0, rank=0, factor_base_size=L + M + 1)

    rows, rhs = [], []
    M_done, Mcur = -1, M
    stage1 = stage2 = 0.0
    anchor = {index[g]: 1}
    with Timer() as total:
        while True:
            # ---------------- Stage 1: sieve every NEW pair (c1 <= c2 <= Mcur, c2 > M_done)
            t0 = time.perf_counter()
            for c1 in range(Mcur + 1):
                lo = max(c1, M_done + 1)
                if lo > Mcur:
                    continue
                A, D, acc = _sieve_segment(c1, lo, Mcur, H, J, primes, logp)
                for i, val in enumerate(acc):
                    c2 = lo + i
                    T = A + c2 * D
                    stats["candidates"] += 1
                    if val < math.log(T) - slack:
                        continue                          # sieve says: not smooth
                    stats["sieve_hits"] += 1
                    exps = smooth_exponents(T, primes, P)  # exact check of the few survivors
                    if exps is None:
                        continue
                    row = {j: e for j, e in enumerate(exps) if e}
                    row[L + c1] = row.get(L + c1, 0) - 1     # - log(H + c1)
                    row[L + c2] = row.get(L + c2, 0) - 1     # - log(H + c2)
                    rows.append(row)
                    rhs.append(0)
                why = limits.exceeded()
                if why:
                    stage1 += time.perf_counter() - t0
                    stats["relations"] = len(rows)
                    raise LinearSieveFailure(f"sieve stopped: {why}", stats)
            M_done = Mcur
            stage1 += time.perf_counter() - t0
            ncols = L + Mcur + 1
            stats.update(relations=len(rows), M_final=Mcur, factor_base_size=ncols)
            # ---------------- Stage 2: sparse elimination mod q, parity lift, verify
            ok = False
            if len(rows) >= 0.9 * ncols:                  # cheap necessary condition
                t0 = time.perf_counter()
                col_count = [0] * ncols
                for row in rows:
                    for c in row:
                        col_count[c] += 1
                col_count[index[g]] += 1
                Xq, rank = solve_mod_prime(rows + [anchor], rhs + [1], ncols, q, col_count)
                stats["rank"] = rank
                ok = all(Xq[i] is not None for i in range(L))
                if ok:
                    logs = []
                    for i, ell in enumerate(primes):
                        lg = lift_to_2q(Xq[i], ell, g, p, q)
                        stats["group_ops"] += 2 * pow_cost(q)
                        if pow(g, lg, p) != ell:
                            raise LinearSieveFailure(f"wrong log for {ell}: inconsistent relations", stats)
                        logs.append(lg)
                stage2 += time.perf_counter() - t0
            if ok:
                break
            if not auto_expand or stats["expansions"] >= max_expand:
                err = LinearSieveFailure(
                    f"insufficient independent relations: {len(rows)} rows (rank {stats['rank']}) "
                    f"for {ncols} unknowns with B={B}, M={Mcur}; increase M/B"
                    + ("" if auto_expand else " or enable auto_expand"), stats)
                err.retryable = True
                raise err
            stats["expansions"] += 1
            Mcur = int(Mcur * growth) + 1
        stats["memory_items"] = len(rows)
        stats["full_column_rank"] = stats["rank"] >= ncols          # rank counts the anchor row too
        stats["factor_base_logs_determined"] = True                 # checked above for every small prime
        # ---------------- Stage 3: individual logarithm over the small primes
        with Timer() as t3:
            x = individual_log(h, g, p, n, q, primes, P, logs, rng, stats, limits)
    stats.update(stage1_sec=stage1, stage2_sec=stage2, stage3_sec=t3.elapsed,
                 time_sec=total.elapsed)
    return x, stats
