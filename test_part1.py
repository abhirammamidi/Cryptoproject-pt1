"""Unit tests for Part 1.   Run:  python test_part1.py   (or: python -m unittest -v test_part1)

* data-file layouts: every supported assignment_data.py layout must normalise to the same instances
* solvers on small generated instances (prime-power orders, composite orders, safe primes)
* the previously failing behaviours (non-invertible rho denominators, independent rho walks, ...)
* a few quick checks on YOUR assignment_data.py (skipped if it has no such instance)
* end to end: run_all.py on a foreign data file (different layout, different numbers) must exit 0
* the primality test against a sieve and against famous strong pseudoprimes
"""
import os
import random
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import bruteforce, bsgs, pollard_rho, pohlig_hellman, index_calculus, linear_sieve
import data_loader
from common import (SolverFailure, element_order, factorize, is_probable_prime, make_instance,
                    primes_up_to, random_safe_prime, safe_prime_generator)


def small_safe_with_g5(bits, rng):
    while True:
        p = random_safe_prime(bits, rng)
        if pow(5, (p - 1) // 2, p) != 1 and pow(5, 2, p) != 1:
            return p


def check(tc, solver, g, h, p, **kw):
    x, st = solver.solve(g, h, p, **kw)
    tc.assertEqual(pow(g, x, p), h)
    tc.assertGreaterEqual(st["time_sec"], 0)
    return x, st


class LoaderTests(unittest.TestCase):
    def setUp(self):
        rng = random.Random(1)
        self.p = random_safe_prime(24, rng)
        self.g = safe_prime_generator(self.p)
        self.h = pow(self.g, 12345, self.p)
        self.tmp = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.tmp.cleanup()

    def load(self, text):
        path = os.path.join(self.tmp.name, "assignment_data.py")
        with open(path, "w") as f:
            f.write(text)
        return data_loader.load_instances(path, quiet=True)[0]

    def expect(self, inst, label="A"):
        self.assertIn(label, inst)
        i = inst[label]
        self.assertEqual((i["p"], i["g"], i["h"]), (self.p, self.g, self.h))
        self.assertEqual(i["order"], self.p - 1)
        self.assertTrue(i["safe_prime"] and i["generator"] and i["h_in_subgroup"])
        self.assertEqual(i["factors"], factorize(self.p - 1))

    def test_dict_of_dicts(self):
        self.expect(self.load(f"INSTANCES={{'A':{{'p':{self.p},'g':{self.g},'h':{self.h}}}}}"))

    def test_tuples_and_other_variable_name(self):
        self.expect(self.load(f"DATA={{'A':({self.p},{self.g},{self.h},{self.p - 1})}}"))

    def test_tuple_in_g_h_p_order(self):
        self.expect(self.load(f"DATA={{'A':({self.g},{self.h},{self.p})}}"))

    def test_variable_per_instance_with_other_field_names(self):
        self.expect(self.load(f"A = {{'P':{self.p},'G':{self.g},'y':{self.h},'N':{self.p - 1}}}"))

    def test_variable_per_field(self):
        self.expect(self.load(f"A_p={self.p}\nA_g={self.g}\nA_h={self.h}\n"))

    def test_list_of_dicts_with_names(self):
        self.expect(self.load(f"X=[{{'name':'Instance A','p':{self.p},'g':{self.g},'h':{self.h}}}]"))

    def test_strings_lowercase_labels(self):
        self.expect(self.load(f"dlp={{'a':{{'prime':'{self.p}','generator':'{self.g}','target':'{self.h}'}}}}"))

    def test_function_hook(self):
        self.expect(self.load(f"def get_instances():\n    return {{'A': ({self.p},{self.g},{self.h})}}\n"))

    def test_supplied_factors_and_wrong_order_are_checked(self):
        inst = self.load(f"INSTANCES={{'A':{{'p':{self.p},'g':{self.g},'h':{self.h},'order':12345}}}}")
        self.assertEqual(inst["A"]["order"], self.p - 1)          # 12345 is not a multiple of ord(g)
        self.assertTrue(inst["A"]["notes"])

    def test_subgroup_order_is_found_even_if_p_minus_1_is_supplied(self):
        q = 1000003
        k = 2
        while not is_probable_prime(k * q + 1):
            k += 1
        p = k * q + 1
        g = next(pow(r, k, p) for r in range(2, 99) if pow(r, k, p) != 1)
        inst = self.load(f"INSTANCES={{'C':{{'p':{p},'g':{g},'h':{pow(g, 77, p)},'order':{p - 1}}}}}")
        self.assertEqual(inst["C"]["order"], q)

    def test_bad_files_give_clear_errors(self):
        with self.assertRaises(data_loader.DataFormatError):
            self.load("X = 1\n")
        with self.assertRaises(data_loader.DataFormatError):
            self.load("INSTANCES = {'A': \n")


class SolverTests(unittest.TestCase):
    def test_bruteforce_and_bsgs_and_rho_and_ph(self):
        rng = random.Random(3)
        for bits in (16, 20, 26):
            g, h, p, n, x = make_instance(bits, rng)
            for solver in (bruteforce, bsgs, pollard_rho, pohlig_hellman):
                got, _ = check(self, solver, g, h, p, order=n)
                self.assertEqual(pow(g, got, p), h)

    def test_edge_exponents(self):
        rng = random.Random(4)
        g, _, p, n, _ = make_instance(20, rng)
        for x in (0, 1, n - 1):
            for solver in (bruteforce, bsgs, pollard_rho, pohlig_hellman):
                check(self, solver, g, pow(g, x, p), p, order=n)

    def test_rho_handles_non_invertible_denominators(self):
        """Even group order (d = 2 candidates) over many seeds: never a wrong answer."""
        p, g = 47, 5
        for h in range(2, 47):
            for seed in range(4):
                check(self, pollard_rho, g, h, p, order=46, seed=seed)

    def test_rho_seeds_give_independent_walks(self):
        rng = random.Random(8)
        g, h, p, n, _ = make_instance(34, rng)
        its = {pollard_rho.solve(g, h, p, order=n, seed=s)[1]["iterations"] for s in range(8)}
        self.assertGreaterEqual(len(its), 7)             # a fixed walk would give (almost) identical counts

    def test_pohlig_hellman_prime_powers_and_inner_counting(self):
        base = 2 * 3**7 * 5**5 * 7**4 * 11**3 * 13**3
        k = 1
        while not is_probable_prime(base * k + 1):
            k += 1
        p = base * k + 1
        facs = factorize(p - 1)
        g = next(c for c in range(2, 99) if element_order(c, p, p - 1, facs) == p - 1)
        h = pow(g, 123456789123456789 % (p - 1), p)
        for inner in ("bsgs", "rho"):
            x, st = pohlig_hellman.solve(g, h, p, order=p - 1, factors=facs, inner=inner)
            self.assertEqual(pow(g, x, p), h)
            self.assertEqual(st["group_ops"], st["outer_group_ops"] + st["inner_group_ops"])
            self.assertGreater(st["inner_group_ops"], 0)

    def test_index_calculus_small_B_does_not_hit_a_trial_cap(self):
        rng = random.Random(9)
        p = small_safe_with_g5(40, rng)
        h = pow(5, rng.randrange(1, p - 1), p)
        for B in (100, 150, 500):                         # B = 100 on a 40-bit prime used to fail
            x, st = check(self, index_calculus, 5, h, p, order=p - 1, B=B, seed=1)
            self.assertGreaterEqual(st["rank"], st["factor_base_size"])

    def test_index_calculus_budget_failure_is_controlled(self):
        rng = random.Random(10)
        p = small_safe_with_g5(40, rng)
        with self.assertRaises(SolverFailure):
            index_calculus.solve(5, pow(5, 99, p), p, order=p - 1, B=100, seed=0, max_trials=10)

    def test_linear_sieve_grows_M_when_parameters_are_too_small(self):
        rng = random.Random(11)
        p = small_safe_with_g5(44, rng)
        h = pow(5, rng.randrange(1, p - 1), p)
        x, st = check(self, linear_sieve, 5, h, p, order=p - 1, B=60, M=60, seed=2)
        self.assertGreater(st["M_final"], st["M_requested"])
        with self.assertRaises(SolverFailure):            # growth switched off: clear message, no hang
            linear_sieve.solve(5, h, p, order=p - 1, B=60, M=60, seed=2, auto_expand=False)

    def test_ic_and_ls_reject_non_generators(self):
        rng = random.Random(12)
        p = random_safe_prime(20, rng)
        with self.assertRaises(ValueError):
            index_calculus.solve(4, pow(4, 7, p), p, order=p - 1, B=50)

    def test_bruteforce_extrapolation_helpers(self):
        rate = bruteforce.calibrate(5, 1000003, 20000)
        self.assertGreater(rate, 0)
        avg_s, worst_s, avg_y, worst_y = bruteforce.extrapolate(rate, 2**64)
        self.assertAlmostEqual(worst_s, 2 * avg_s, delta=avg_s * 1e-9)


class OwnDataTests(unittest.TestCase):
    """Quick checks on the real assignment_data.py (any instance that is missing is skipped)."""
    @classmethod
    def setUpClass(cls):
        try:
            cls.inst, _ = data_loader.load_instances(quiet=True)
        except data_loader.DataFormatError as exc:
            raise unittest.SkipTest(str(exc)) from exc

    def need(self, label):
        if label not in self.inst:
            self.skipTest(f"instance {label} not in the data file")
        return self.inst[label]

    def run_solver(self, solver, label, **kw):
        i = self.need(label)
        check(self, solver, i["g"], i["h"], i["p"], order=i["order"], **kw)

    def test_A_bruteforce_bsgs(self):
        self.run_solver(bruteforce, "A")
        self.run_solver(bsgs, "A")

    def test_D_E_pohlig_hellman(self):
        for lab in ("D", "E"):
            self.run_solver(pohlig_hellman, lab, factors=self.need(lab)["factors"])

    def test_F_index_calculus(self):
        self.run_solver(index_calculus, "F", B=500, seed=1)

    def test_G_linear_sieve(self):
        self.run_solver(linear_sieve, "G", B=700, M=900, seed=1)


def small_dataset(rng):
    """A complete A-H(+C2) set of SMALL instances with the right structure (seconds to solve)."""
    def rand_prime(bits):
        while True:
            c = rng.getrandbits(bits) | (1 << (bits - 1)) | 1
            if is_probable_prime(c):
                return c

    def with_order(p, n):                          # element of exact order n
        facs = factorize(n)
        while True:
            g = pow(rng.randrange(2, p - 1), (p - 1) // n, p)
            if g != 1 and element_order(g, p, n, facs) == n:
                return g

    def prime_with_factor(q):                      # p = k*q + 1
        k = 2
        while not is_probable_prime(k * q + 1):
            k += 2
        return k * q + 1

    data = {}
    for lab, bits in (("A", 16), ("B", 20)):
        g, h, p, n, _ = make_instance(bits, rng)
        data[lab] = (p, g, h)
    q = rand_prime(22)
    p = prime_with_factor(q)
    data["C"] = (p, with_order(p, q), None)
    q2 = rand_prime(12) * rand_prime(12)
    p = prime_with_factor(q2)
    data["C2"] = (p, with_order(p, q2), None)
    small = primes_up_to(200)[1:]
    while True:                                    # D: smooth order
        n = 2
        while n.bit_length() < 44:
            n *= rng.choice(small)
        if is_probable_prime(n + 1):
            break
    data["D"] = (n + 1, with_order(n + 1, n), None)
    base, k = 2 * 3**4 * 5**3 * 7**2, 1            # E: prime powers
    while not is_probable_prime(base * k + 1):
        k += 1
    data["E"] = (base * k + 1, with_order(base * k + 1, base * k), None)
    for lab, bits in (("F", 28), ("G", 32), ("H", 34)):
        p = small_safe_with_g5(bits, rng)
        data[lab] = (p, 5, None)
    out = {}
    for lab, (p, g, h) in data.items():
        if h is None:
            n = element_order(g, p, p - 1, factorize(p - 1))
            h = pow(g, rng.randrange(1, n), p)
        out[lab] = (p, g, h)
    return out


class EndToEndTests(unittest.TestCase):
    """run_all.py as the examiner would run it, on a data file it has never seen."""
    def run_all(self, text, *extra):
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "assignment_data.py")
            with open(path, "w") as f:
                f.write(text)
            return subprocess.run([sys.executable, os.path.join(HERE, "run_all.py"), "--data", path,
                                   "--out", os.path.join(d, "results"), *extra],
                                  capture_output=True, text=True, timeout=900, cwd=d)

    def test_foreign_data_file_runs_cleanly(self):
        data = small_dataset(random.Random(21))
        r = self.run_all("PROBLEM_SET = " + repr(data) + "\n", "--required-only")   # tuples, another name
        self.assertEqual(r.returncode, 0, r.stdout[-3000:] + r.stderr[-3000:])
        self.assertIn("TABLE 2", r.stdout)
        self.assertIn("0 required failures", r.stdout)
        self.assertIn("12 solved and verified", r.stdout)                     # all 12 required solves

    def test_unreadable_data_file_gives_message_not_traceback(self):
        r = self.run_all("X = 1\n")
        self.assertEqual(r.returncode, 2)
        self.assertIn("ERROR", r.stdout)
        self.assertNotIn("Traceback", r.stdout + r.stderr)


class PrimalityTests(unittest.TestCase):
    def test_against_a_sieve(self):
        primes = set(primes_up_to(30000))
        for n in range(30000):
            self.assertEqual(is_probable_prime(n), n in primes, n)

    def test_strong_pseudoprimes_are_rejected(self):
        # composites that fool Miller-Rabin for many fixed bases (the last one fools bases 2..37)
        for n in (2047, 3215031751, 2152302898747, 3474749660383, 341550071728321,
                  3825123056546413051, 318665857834031151167461, 3317044064679887385961981):
            self.assertFalse(is_probable_prime(n), n)
        for n in (2**61 - 1, 2**89 - 1, 2**127 - 1, 1000000007):
            self.assertTrue(is_probable_prime(n), n)


if __name__ == "__main__":
    unittest.main(verbosity=2)
