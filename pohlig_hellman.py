"""Task 1.3 - Pohlig-Hellman with prime-power lifting and CRT."""

import time
from common import factorize, crt
import bsgs


def _dlog_order_q(g, h, p, q):
    return bsgs.solve(g, h, p, order=q)[0]


def _prime_power_log(g, h, p, n, q, e, stats):
    qe = q ** e

    # Project g and h into the subgroup of order q^e.
    gi = pow(g, n // qe, p)
    hi = pow(h, n // qe, p)

    # gamma has order q.
    gamma = pow(gi, q ** (e - 1), p)
    gi_inv = pow(gi, qe - 1, p)

    x = 0

    for k in range(e):
        # Remove the digits already recovered.
        residual = hi * pow(gi_inv, x, p) % p

        # Project the remaining unknown digit into the order-q subgroup.
        target = pow(residual, q ** (e - 1 - k), p)

        if target == 1:
            digit = 0
        else:
            digit = _dlog_order_q(gamma, target, p, q)
            if digit is None:
                raise ValueError("inner DLP failed during digit lifting")

        x += digit * (q ** k)
        stats["largest_subgroup"] = max(stats["largest_subgroup"], q)

    return x


def solve(g, h, p, order=None, factors=None, **kwargs):
    n = p - 1 if order is None else order
    if factors is None:
        factors = factorize(n)

    stats = {
        "group_ops": 0, "time_sec": 0.0, "memory_items": 1,
        "largest_subgroup": 1, "factors": dict(sorted(factors.items()))
    }

    start = time.perf_counter()
    residues = []
    moduli = []

    for q, e in sorted(factors.items()):
        r = _prime_power_log(g, h, p, n, q, e, stats)
        residues.append(r)
        moduli.append(q ** e)

    x, combined = crt(residues, moduli)

    stats["time_sec"] = time.perf_counter() - start

    if pow(g, x, p) != h:
        raise ValueError("Pohlig-Hellman verification failed")

    return x, stats


if __name__ == "__main__":
    x, st = solve(5, 8, 23, 22, factors={2: 1, 11: 1})
    print("x =", x)
    print("verification =", pow(5, x, 23) == 8)
    print(st)
