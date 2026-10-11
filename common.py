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


class SolverFailure(RuntimeError):
    """A solver gave up in a controlled way (budget, rank, parameters).
    `.stats` carries diagnostics so a failed run is still a documented result."""
    def __init__(self, message, stats=None):
        super().__init__(message)
        self.stats = stats or {}


class Limits:
    """Wall-clock and trial budgets. exceeded(trials) returns a reason or None."""
    def __init__(self, max_trials=None, max_seconds=None):
        self.max_trials, self.max_seconds = max_trials, max_seconds
        self.t0 = time.perf_counter()

    def exceeded(self, trials=0):
        if self.max_trials is not None and trials > self.max_trials:
            return f"trial budget {self.max_trials:,} exhausted"
        if self.max_seconds is not None and time.perf_counter() - self.t0 > self.max_seconds:
            return f"time budget {self.max_seconds:g}s exhausted"
        return None


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


# Miller-Rabin with the first 13 prime bases is PROVEN correct for every n below this bound
# (Sorenson & Webster, 2015). Above it, 20 extra pseudo-random bases are added.
_MR_BASES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41)
_MR_PROVEN_BELOW = 3317044064679887385961981


def is_probable_prime(n, extra_rounds=20):
    """Miller-Rabin. Deterministic (a proof, not a probability) for n < 3.3 * 10^24 (~81 bits);
    for larger n the error probability is below 4^-(13 + extra_rounds)."""
    if n < 2:
        return False
    for q in _SMALL[:25]:
        if n % q == 0:
            return n == q
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    bases = list(_MR_BASES)
    if n >= _MR_PROVEN_BELOW:
        rng = random.Random(n)                    # reproducible extra bases
        bases += [rng.randrange(2, n - 1) for _ in range(extra_rounds)]
    for a in bases:
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


def _pollard_brent(n, rng, deadline):
    """Brent's variant of Pollard's rho FACTORING method: one non-trivial factor of composite n."""
    if n % 2 == 0:
        return 2
    while True:
        y, c, m = rng.randrange(1, n), rng.randrange(1, n), 128
        g = r = q = 1
        x = ys = y
        while g == 1:
            x = y
            for _ in range(r):
                y = (y * y + c) % n
            k = 0
            while k < r and g == 1:
                ys = y
                for _ in range(min(m, r - k)):
                    y = (y * y + c) % n
                    q = q * abs(x - y) % n
                g = math.gcd(q, n)
                k += m
            r *= 2
            if time.perf_counter() > deadline:
                raise ValueError("factorisation time budget exhausted")
        if g == n:
            g = 1
            while g == 1:
                ys = (ys * ys + c) % n
                g = math.gcd(abs(x - ys), n)
        if g != n:
            return g


def factorize(n, max_seconds=30):
    """Return {prime: exponent}: trial division by small primes, then Miller-Rabin and
    Pollard-Brent rho on what is left. Raises ValueError if it cannot finish in time."""
    if n < 1:
        raise ValueError("cannot factor a non-positive number")
    factors = {}
    for q in _SMALL:
        while n % q == 0:
            factors[q] = factors.get(q, 0) + 1
            n //= q
    deadline = time.perf_counter() + max_seconds
    rng = random.Random(n)
    stack = [n] if n > 1 else []
    while stack:
        m = stack.pop()
        if m == 1:
            continue
        if is_probable_prime(m):
            factors[m] = factors.get(m, 0) + 1
            continue
        d = _pollard_brent(m, rng, deadline)
        stack += [d, m // d]
    return dict(sorted(factors.items()))


def available_ram_bytes(default=4 * 2**30):
    """Best-effort physical RAM size (Linux/macOS via sysconf, Windows via ctypes)."""
    try:
        import os
        return os.sysconf("SC_PHYS_PAGES") * os.sysconf("SC_PAGE_SIZE")
    except (ValueError, OSError, AttributeError):
        pass
    try:
        import ctypes

        class _MS(ctypes.Structure):
            _fields_ = [("l", ctypes.c_ulong), ("load", ctypes.c_ulong), ("total", ctypes.c_ulonglong),
                        ("avail", ctypes.c_ulonglong), ("a", ctypes.c_ulonglong), ("b", ctypes.c_ulonglong),
                        ("c", ctypes.c_ulonglong), ("d", ctypes.c_ulonglong), ("e", ctypes.c_ulonglong)]
        ms = _MS()
        ms.l = ctypes.sizeof(_MS)
        ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(ms))
        return int(ms.total)
    except Exception:
        return default


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
