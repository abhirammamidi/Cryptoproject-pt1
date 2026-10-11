"""Task 1.3 - Pohlig-Hellman with prime-power digit lifting and CRT."""
from common import new_stats, Timer, pow_cost, factorize, crt
import bsgs
import pollard_rho

INNER_SOLVERS = {"bsgs": bsgs.solve, "rho": pollard_rho.solve}


def _prime_power_log(g, h, p, n, q, e, inner, stats):
    """Return x mod q^e, where g has order n and q^e exactly divides n."""
    qe = q ** e
    gi = pow(g, n // qe, p)                 # order q^e
    hi = pow(h, n // qe, p)
    gamma = pow(gi, q ** (e - 1), p)        # order q (the digit base)
    stats["group_ops"] += pow_cost(n // qe) * 2 + pow_cost(q ** (e - 1))
    gi_inv = pow(gi, qe - 1, p)             # gi^(-1), since gi^(q^e) = 1
    stats["group_ops"] += pow_cost(qe - 1)

    x = 0
    for k in range(e):
        # remove the digits already known, then push into the order-q subgroup
        stripped = hi * pow(gi_inv, x, p) % p
        hk = pow(stripped, q ** (e - 1 - k), p)
        stats["group_ops"] += pow_cost(x) + 1 + pow_cost(q ** (e - 1 - k))
        if hk == 1:
            dk = 0                          # digit is zero, no DLP needed
        else:
            dk, sub = inner(gamma, hk, p, order=q, seed=k)
            stats["group_ops"] += sub["group_ops"]
            stats["memory_items"] = max(stats["memory_items"], sub["memory_items"])
        x += dk * q ** k
        stats["largest_subgroup"] = max(stats["largest_subgroup"], q)
    return x


def solve(g, h, p, order=None, factors=None, inner="bsgs", **kwargs):
    n = (p - 1) if order is None else order
    stats = new_stats()
    stats["largest_subgroup"] = 1
    inner_solve = INNER_SOLVERS[inner]
    with Timer() as t:
        if factors is None:
            factors = factorize(n)
        residues, moduli = [], []
        for q, e in sorted(factors.items()):
            residues.append(_prime_power_log(g, h, p, n, q, e, inner_solve, stats))
            moduli.append(q ** e)
        x, _ = crt(residues, moduli)
    stats["time_sec"] = t.elapsed
    stats["factors"] = dict(sorted(factors.items()))
    if pow(g, x, p) != h:
        raise ValueError("Pohlig-Hellman failed: check order/factorisation")
    return x, stats
