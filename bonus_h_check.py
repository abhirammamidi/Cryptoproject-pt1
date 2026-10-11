"""Verify the supplied H bonus test secret.

IMPORTANT: This does NOT solve the discrete logarithm. It reads the known
verification exponent from assignment_data.TEST_SECRETS and checks it against
H's public parameters. Do not describe this as a Pollard-rho/algorithmic solve.
"""
import assignment_data


def main():
    inst = assignment_data.INSTANCES["H"]
    x = assignment_data.TEST_SECRETS["H"]
    lhs = pow(inst["g"], x, inst["p"])
    ok = lhs == inst["h"]
    print("H supplied test-secret verification")
    print(f"p = {inst['p']}")
    print(f"g = {inst['g']}")
    print(f"h = {inst['h']}")
    print(f"x = {x}")
    print(f"g^x mod p = {lhs}")
    print(f"verified = {ok}")
    if not ok:
        raise SystemExit("H test secret did not verify")
    print("NOTE: x was read from TEST_SECRETS; no DLP algorithm recovered it.")


if __name__ == "__main__":
    main()
