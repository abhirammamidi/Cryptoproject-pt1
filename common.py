"""Shared helpers for Part 1 reference implementations.

No discrete-log library is used here.  Core algorithms use only the Python
standard library.
"""
import math
import random
import time


def new_stats():
    return {"group_ops": 0, "time_sec": 0.0, "memory_items": 0}


class Timer:
    def __enter__(self):
        self.start = time.perf_counter()
        return self

    def __exit__(self, *_):
        self.elapsed = time.perf_counter() - self.start
        return False


def pow_cost(e: int) -> int:
    """Approximate modular multiplications for square-and-multiply."""
    if e <= 1:
        return 0
    return (e.bit_length() - 1) + (e.bit_count() - 1)


def egcd(a, b):
    old_r, r = a, b
    old_s, s = 1, 0
    old_t, t = 0, 1
    while r:
        q = old_r // r
        old_r, r = r, old_r - q * r
        old_s, s = s, old_s - q * s
        old_t, t = t, old_t - q * t
    return old_r, old_s, old_t


def modinv(a, m):
    g, x, _ = egcd(a % m, m)
    if g != 1:
        raise ValueError("inverse does not exist")
    return x % m


def solve_linear_congruence(a, b, m):
    """All x in [0,m) satisfying a*x == b (mod m)."""
    a %= m
    b %= m
    d = math.gcd(a, m)
    if b % d:
        return []
    a1, b1, m1 = a // d, b // d, m // d
    x0 = (b1 * modinv(a1, m1)) % m1
    return [x0 + k * m1 for k in range(d)]


def crt(residues, moduli):
    x, M = 0, 1
    for r, m in zip(residues, moduli):
        t = ((r - x) * modinv(M % m, m)) % m
        x += M * t
        M *= m
        x %= M
    return x, M


def primes_up_to(limit):
    if limit < 2:
        return []
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[0:2] = b"\x00\x00"
    for p in range(2, math.isqrt(limit) + 1):
        if sieve[p]:
            sieve[p*p:limit+1:p] = b"\x00" * (((limit - p*p)//p) + 1)
    return [i for i, ok in enumerate(sieve) if ok]


_SMALL_PRIMES = primes_up_to(1000)


def is_probable_prime(n, rounds=16):
    if n < 2:
        return False
    small = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
    for p in small:
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        s += 1
        d //= 2
    # Deterministic bases for common 64-bit-sized test data.
    bases = (2, 325, 9375, 28178, 450775, 9780504, 1795265022)
    for a in bases[:max(1, min(len(bases), rounds))]:
        a %= n
        if a in (0, 1):
            continue
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


def factorize(n):
    """Trial division; intended for assignment orders whose factorization is supplied."""
    out = {}
    for p in _SMALL_PRIMES:
        while n % p == 0:
            out[p] = out.get(p, 0) + 1
            n //= p
        if p * p > n:
            break
    if n > 1:
        if not is_probable_prime(n):
            raise ValueError("remaining cofactor is not prime; supply factors explicitly")
        out[n] = out.get(n, 0) + 1
    return out


def safe_prime_generator(p):
    """Find a generator of F_p* when p=2q+1 and q is prime."""
    q = (p - 1) // 2
    for g in range(2, min(p, 10000)):
        if pow(g, 2, p) != 1 and pow(g, q, p) != 1:
            return g
    raise ValueError("generator search failed")


def random_safe_prime(bits, rng):
    while True:
        q = rng.getrandbits(bits - 1) | (1 << (bits - 2)) | 1
        if is_probable_prime(q) and is_probable_prime(2*q + 1):
            return 2*q + 1
