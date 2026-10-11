"""Run Part 1 against assignment_data.py.

Expected assignment_data.py interface:
INSTANCES = {
    "A": {"p": ..., "g": ..., "h": ..., "order": ...},
    ...
}
Optional: factors, B_ic, B_ls, M_ls.
"""

import argparse
import importlib
import math

import bruteforce
import bsgs
import pollard_rho
import pohlig_hellman
import index_calculus
import linear_sieve


def load_instances():
    data = importlib.import_module("assignment_data")
    if not hasattr(data, "INSTANCES"):
        raise ValueError("assignment_data.py must define INSTANCES")
    return data.INSTANCES


def verify(inst, x):
    return x is not None and pow(inst["g"], x, inst["p"]) == inst["h"]


SOLVERS = {
    "bruteforce": lambda i: bruteforce.solve(
        i["g"], i["h"], i["p"], i.get("order")
    ),
    "bsgs": lambda i: bsgs.solve(
        i["g"], i["h"], i["p"], i.get("order")
    ),
    "rho": lambda i: pollard_rho.solve(
        i["g"], i["h"], i["p"], i.get("order"), seed=0
    ),
    "pohlig_hellman": lambda i: pohlig_hellman.solve(
        i["g"], i["h"], i["p"], i.get("order"),
        factors=i.get("factors")
    ),
    "index_calculus": lambda i: index_calculus.solve(
        i["g"], i["h"], i["p"], i.get("order"),
        B=i.get("B_ic", 500), seed=0
    ),
    "linear_sieve": lambda i: linear_sieve.solve(
        i["g"], i["h"], i["p"], i.get("order"),
        B=i.get("B_ls", 700), M=i.get("M_ls", 900), seed=0
    ),
}


def estimate(algorithm, inst):
    n = inst.get("order", inst["p"]-1)
    if algorithm == "bruteforce":
        return n/2
    if algorithm in ("bsgs", "rho"):
        return 2*math.isqrt(n)
    if algorithm == "pohlig_hellman":
        factors = inst.get("factors")
        if factors:
            return 2 * sum(math.isqrt(q) * e for q, e in factors.items())
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default="bruteforce,bsgs,rho,pohlig_hellman,index_calculus,linear_sieve")
    ap.add_argument("--budget-ops", type=float, default=2e8)
    args = ap.parse_args()

    instances = load_instances()

    for name in args.only.split(","):
        if name not in SOLVERS:
            print("Unknown solver:", name)
            continue

        print(f"\n=== {name} ===")
        for label in "ABCDEFGH":
            if label not in instances:
                continue

            inst = instances[label]
            est = estimate(name, inst)
            if est is not None and est > args.budget_ops:
                print(f"{label}: SKIPPED (estimated {est:.3e} operations)")
                continue

            try:
                x, stats = SOLVERS[name](inst)
                ok = verify(inst, x)
                print(f"{label}: x={x}  verified={ok}  stats={stats}")
            except Exception as exc:
                print(f"{label}: FAILED: {type(exc).__name__}: {exc}")


if __name__ == "__main__":
    main()
