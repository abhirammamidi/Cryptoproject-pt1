"""Shared number-theory and measurement helpers for Part 1.

Everything here is standard-library only. None of these functions solves a
discrete logarithm; they are the auxiliary tools the six solvers rely on.
"""
import math
import random
import time


# ---------------------------------------------------------------- counting
def pow_cost(e):
    """Group multiplications used by left-to-right square-and-multiply
    to compute a**e: (bit_length - 1) squarings + (popcount - 1) multiplies."""
    if e <= 1:
        return 0
    return (e.bit_length() - 1) + (bin(e).count("1") - 1)


def new_stats():
    return {"group_ops": 0, "time_sec": 0.0, "memory_items": 0}


class Timer:
    """Usage:  with Timer() as t: ...   then read t.elapsed (seconds)."""
    def __enter__(self):
        self._start = time.perf_counter()
        return self

    def __exit__(self, *exc):
        self.elapsed = time.perf_counter() - self._start
        return False


# ---------------------------------------------------------------- arithmetic
def egcd(a, b):
    """Return (g, s, t) with a*s + b*t == g == gcd(a, b)."""
    old_r, r = a, b
    old_s, s = 1, 0
    old_t, t = 0, 1
    while r != 0:
        q = old_r // r
        old_r, r = r, old_r - q * r
        old_s, s = s, old_s - q * s
        old_t, t = t, old_t - q * t
    return old_r, old_s, old_t


def modinv(a, m):
    g, s, _ = egcd(a % m, m)
    if g != 1:
        raise ValueError(f"{a} has no inverse mod {m} (gcd = {g})")
    return s % m


def solve_linear_congruence(a, b, m):
    """All x in [0, m) with a*x == b (mod m). May return an empty list."""
    a, b = a % m, b % m
    d, _, _ = egcd(a, m)
    if b % d != 0:
        return []
    a1, b1, m1 = a // d, b // d, m // d
    x0 = (b1 * modinv(a1, m1)) % m1
    return [x0 + k * m1 for k in range(d)]


def crt(residues, moduli):
    """Combine x = r_i (mod m_i) for pairwise-coprime m_i. Returns (x, M)."""
    X, M = 0, 1
    for r, m in zip(residues, moduli):
        t = ((r - X) * modinv(M % m, m)) % m
        X, M = X + M * t, M * m
    return X % M, M


# ---------------------------------------------------------------- primes
def primes_up_to(B):
    """Sieve of Eratosthenes: all primes <= B (B itself included if prime)."""
    if B < 2:
        return []
    sieve = bytearray([1]) * (B + 1)
    sieve[0] = sieve[1] = 0
    for i in range(2, math.isqrt(B) + 1):
        if sieve[i]:
            sieve[i * i::i] = bytearray(len(range(i * i, B + 1, i)))
    return [i for i in range(B + 1) if sieve[i]]


_SMALL = primes_up_to(1000)


def is_probable_prime(n, rounds=32):
    """Miller-Rabin primality test."""
    if n < 2:
        return False
    for q in _SMALL[:25]:
        if n % q == 0:
            return n == q
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    rng = random.Random(n)            # deterministic per n
    for _ in range(rounds):
        a = rng.randrange(2, n - 1)
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def factorize(n, trial_limit=10**6):
    """Return {prime: exponent}. Trial division up to trial_limit; a leftover
    cofactor is accepted only if Miller-Rabin says it is prime."""
    factors = {}
    for q in _SMALL:
        while n % q == 0:
            factors[q] = factors.get(q, 0) + 1
            n //= q
    q = _SMALL[-1] + 2
    while n > 1 and q * q <= n and q <= trial_limit:
        while n % q == 0:
            factors[q] = factors.get(q, 0) + 1
            n //= q
        q += 2
    if n > 1:
        if not is_probable_prime(n):
            raise ValueError(f"could not fully factor; composite cofactor {n}")
        factors[n] = factors.get(n, 0) + 1
    return factors


def element_order(g, p, n=None, factors=None):
    """Exact order of g mod p, given a multiple n of it (default p - 1)."""
    n = p - 1 if n is None else n
    factors = factorize(n) if factors is None else factors
    order = n
    for q, e in factors.items():
        for _ in range(e):
            if pow(g, order // q, p) == 1:
                order //= q
            else:
                break
    return order


# ---------------------------------------------------------------- test data
def random_safe_prime(bits, rng):
    """p = 2q + 1 with p and q prime and p exactly `bits` bits long."""
    while True:
        q = rng.getrandbits(bits - 1) | (1 << (bits - 2)) | 1
        if is_probable_prime(q) and is_probable_prime(2 * q + 1):
            return 2 * q + 1


def safe_prime_generator(p):
    """Smallest g of order p - 1 for a safe prime p = 2q + 1."""
    q = (p - 1) // 2
    for g in range(2, p):
        if pow(g, 2, p) != 1 and pow(g, q, p) != 1:
            return g


def make_instance(bits, rng):
    """A random toy DLP (g, h, p, n, x) on a safe prime - for self-testing
    and experiments only, never a substitute for the official instances."""
    p = random_safe_prime(bits, rng)
    g = safe_prime_generator(p)
    n = p - 1
    x = rng.randrange(1, n)
    return g, pow(g, x, p), p, n, x
