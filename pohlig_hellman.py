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
    stats["outer_group_ops"] += pow_cost(n // qe) * 2 + pow_cost(q ** (e - 1))
    gi_inv = pow(gi, qe - 1, p)             # gi^(-1), since gi^(q^e) = 1
    stats["outer_group_ops"] += pow_cost(qe - 1)

    x, digits = 0, []
    for k in range(e):
        # remove the digits already known, then push into the order-q subgroup
        stripped = hi * pow(gi_inv, x, p) % p
        hk = pow(stripped, q ** (e - 1 - k), p)
        stats["outer_group_ops"] += pow_cost(x) + 1 + pow_cost(q ** (e - 1 - k))
        if hk == 1:
            dk = 0                          # digit is zero, no DLP needed
        else:
            dk, sub = inner(gamma, hk, p, order=q, seed=k)
            stats["inner_group_ops"] += sub["group_ops"]          # inner solver work IS counted
            stats["inner_peak_table"] = max(stats["inner_peak_table"], sub["memory_items"])
            stats["inner_calls"] += 1
        digits.append(dk)
        x += dk * q ** k
        stats["largest_subgroup"] = max(stats["largest_subgroup"], q)
    stats["subproblems"].append({"q": q, "e": e, "digits": digits, "x_mod_qe": x})
    return x


def solve(g, h, p, order=None, factors=None, inner="bsgs", **kwargs):
    """Pohlig-Hellman: one small DLP per prime power q^e of n (digit lifting), then CRT.

    Args:
        g, h, p: the problem g^x = h (mod p), p prime.
        order:   order n of g (default p - 1); the answer is returned modulo n.
        factors: {q: e} with n = product of q^e (computed with common.factorize if omitted).
        inner:   "bsgs" (default) or "rho" for the order-q sub-problems.
    Returns:
        (x, stats) with stats = {group_ops, time_sec, memory_items, outer_group_ops,
        inner_group_ops, inner_peak_table, inner_calls, largest_subgroup, subproblems, factors}.
        group_ops = outer exponentiations + every inner solver's group_ops;
        memory_items = the largest inner BSGS table.
    Raises:
        ValueError if the factorisation does not match the order (the CRT result fails to verify).
    """
    n = (p - 1) if order is None else order
    stats = new_stats()
    stats.update(largest_subgroup=1, outer_group_ops=0, inner_group_ops=0,
                 inner_peak_table=0, inner_calls=0, subproblems=[])
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
    stats["group_ops"] = stats["outer_group_ops"] + stats["inner_group_ops"]
    stats["memory_items"] = stats["inner_peak_table"]
    stats["factors"] = dict(sorted(factors.items()))
    if pow(g, x, p) != h:
        raise ValueError("Pohlig-Hellman failed: check order/factorisation")
    return x, stats
