import unittest
import assignment_data
import bruteforce, bsgs, pollard_rho, pohlig_hellman, index_calculus, linear_sieve

class SolverTests(unittest.TestCase):
    def check(self, solver, label, **kwargs):
        inst = assignment_data.INSTANCES[label]
        x, stats = solver.solve(inst["g"], inst["h"], inst["p"],
                                order=inst["order"], **kwargs)
        self.assertEqual(pow(inst["g"], x, inst["p"]), inst["h"])
        self.assertGreaterEqual(stats["time_sec"], 0)
        return x

    def test_bruteforce_A(self):
        self.check(bruteforce, "A")

    def test_bsgs_A(self):
        self.check(bsgs, "A")

    def test_rho_C2(self):
        self.check(pollard_rho, "C2", seed=12345, max_restarts=100)

    def test_ph_E(self):
        i = assignment_data.INSTANCES["E"]
        self.check(pohlig_hellman, "E", factors=i["factors"])

    def test_ic_F(self):
        self.check(index_calculus, "F", B=500, seed=123)

    def test_linear_sieve_G(self):
        x = self.check(linear_sieve, "G", B=700, M=900, seed=20261011)
        self.assertEqual(x, assignment_data.TEST_SECRETS["G"])

if __name__ == "__main__":
    unittest.main()
