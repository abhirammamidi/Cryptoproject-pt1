"""Run every algorithm assigned to each DLP instance and print a results table.

The runner follows assignment_data.INSTANCES[*]["tasks"] so it does not launch
computationally infeasible or mathematically incompatible algorithms on every
input. Use --include-extra only if you intentionally want those extra trials.
"""
from __future__ import annotations
import argparse
import csv
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import assignment_data
import bruteforce, bsgs, index_calculus, linear_sieve
import pohlig_hellman, pollard_rho

RESULTS = HERE / "experiments" / "results"
SOLVERS = {
    "bruteforce": bruteforce.solve,
    "bsgs": bsgs.solve,
    "rho": pollard_rho.solve,
    "pohlig_hellman": pohlig_hellman.solve,
    "index_calculus": index_calculus.solve,
    "linear_sieve": linear_sieve.solve,
}
DISPLAY = {
    "bruteforce": "Brute force",
    "bsgs": "Baby-step/Giant-step",
    "rho": "Pollard rho",
    "pohlig_hellman": "Pohlig-Hellman",
    "index_calculus": "Index calculus",
    "linear_sieve": "Linear sieve",
    "bonus": "Exponent validation",
}

def solver_kwargs(label, inst, algorithm, seed):
    kw = {"order": inst["order"]}
    if algorithm == "pohlig_hellman":
        kw["factors"] = inst["factors"]
    elif algorithm == "rho":
        kw.update(seed=seed, max_restarts=100)
    elif algorithm == "index_calculus":
        kw.update(B={"F": 500, "G": 2000}.get(label, inst.get("B_ic", 500)),
                  seed=seed)
    elif algorithm == "linear_sieve":
        kw.update(B=700, M=900, seed=seed)
    return kw

def counted_modexp(base, exponent, modulus):
    """Square-and-multiply verification with a count of modular multiplications."""
    result = 1
    base %= modulus
    group_ops = 0
    while exponent > 0:
        if exponent & 1:
            result = (result * base) % modulus
            group_ops += 1
        exponent >>= 1
        if exponent:
            base = (base * base) % modulus
            group_ops += 1
    return result, group_ops


def run_one(label, inst, algorithm, seed):
    x, stats = SOLVERS[algorithm](
        inst["g"], inst["h"], inst["p"],
        **solver_kwargs(label, inst, algorithm, seed)
    )
    if pow(inst["g"], x, inst["p"]) != inst["h"]:
        raise AssertionError("independent verification failed")
    stats = dict(stats)
    stats["verified"] = True
    return x, stats

def fmt(value, digits=6):
    if value in ("", None):
        return "-"
    try:
        return f"{float(value):.{digits}f}"
    except (TypeError, ValueError):
        return str(value)

def print_table(rows):
    headers = ["Instance", "Algorithm", "Status", "Recovered x", "Time (s)", "Group ops", "Memory"]
    data = []
    for r in rows:
        x = str(r.get("x", "-"))
        data.append([
            r["instance"], r["algorithm"], r["status"], x,
            fmt(r.get("time_sec")), str(r.get("group_ops", "-")),
            str(r.get("memory_items", "-"))
        ])
    widths = [len(h) for h in headers]
    for row in data:
        widths = [max(w, len(str(v))) for w, v in zip(widths, row)]
    line = "+-" + "-+-".join("-" * w for w in widths) + "-+"
    print("\n" + line)
    print("| " + " | ".join(h.ljust(w) for h, w in zip(headers, widths)) + " |")
    print(line)
    for row in data:
        print("| " + " | ".join(str(v).ljust(w) for v, w in zip(row, widths)) + " |")
    print(line)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=20261011)
    parser.add_argument("--only", default=None,
                        help="comma-separated algorithm names: " + ", ".join(SOLVERS))
    parser.add_argument("--include-extra", action="store_true",
                        help="try every selected solver on every instance; may be very slow or fail")
    args = parser.parse_args()
    selected = set(args.only.split(",")) if args.only else set(SOLVERS)
    unknown = selected - set(SOLVERS)
    if unknown:
        parser.error("unknown algorithm(s): " + ", ".join(sorted(unknown)))

    rows = []
    failures = 0
    print("Discrete Logarithm Assignment — Part 1")
    print("Running assigned algorithms; each successful x is independently verified.")
    for label, inst in assignment_data.INSTANCES.items():
        tasks = list(SOLVERS) if args.include_extra else list(inst.get("tasks", ()))
        for algorithm in tasks:
            if algorithm == "bonus":
                # This is explicitly not counted as a DLP solver result.
                x_known = assignment_data.TEST_SECRETS.get(label)
                started = time.perf_counter()
                if x_known is not None:
                    check_value, verification_ops = counted_modexp(
                        inst["g"], x_known, inst["p"])
                    verified = check_value == inst["h"]
                else:
                    verification_ops = 0
                    verified = False
                verification_time = time.perf_counter() - started
                rows.append({
                    "instance": label, "algorithm": DISPLAY["bonus"],
                    "status": "PASS" if verified else "FAIL",
                    # Supplied verification secret, not a recovered DLP result.
                    "x": x_known if verified else "-",
                    # Count modular multiplications used by square-and-multiply.
                    "time_sec": verification_time,
                    "group_ops": verification_ops,
                    # One supplied secret value is retained; this is a simple item count,
                    # not a byte-accurate or peak-memory measurement.
                    "memory_items": 1 if x_known is not None else 0
                })
                print(f"[{label}] Exponent validation "
                      f"{'PASS' if verified else 'FAIL'}; "
                      f"{verification_time:.6f}s, {verification_ops} modular multiplications.")
                if not verified:
                    failures += 1
                continue
            if algorithm not in selected:
                continue
            seed = args.seed + sum((i + 1) * ord(c)
                                   for i, c in enumerate(label + algorithm))
            print(f"[{label}] {DISPLAY[algorithm]} ...", flush=True)
            try:
                x, stats = run_one(label, inst, algorithm, seed)
                rows.append({
                    "instance": label, "algorithm": DISPLAY[algorithm],
                    "status": "PASS", "x": x,
                    "time_sec": stats.get("time_sec", 0),
                    "group_ops": stats.get("group_ops", "-"),
                    "memory_items": stats.get("memory_items", "-"),
                })
                print(f"    recovered and verified (time {stats.get('time_sec', 0):.4f}s)")
            except Exception as exc:
                failures += 1
                rows.append({
                    "instance": label, "algorithm": DISPLAY[algorithm],
                    "status": "FAIL: " + type(exc).__name__,
                    "x": "-", "time_sec": "-", "group_ops": "-",
                    "memory_items": "-"
                })
                print(f"    FAILED: {type(exc).__name__}: {exc}")

    print_table(rows)
    RESULTS.mkdir(parents=True, exist_ok=True)
    csv_path = RESULTS / "results.csv"
    fields = ["instance", "algorithm", "status", "x", "time_sec", "group_ops", "memory_items"]
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    with (RESULTS / "summary.md").open("w", encoding="utf-8") as f:
        f.write("# Part 1 results\n\n")
        f.write("| Instance | Algorithm | Status | Recovered x | Time (s) | Group ops | Memory |\n")
        f.write("|---|---|---|---:|---:|---:|---:|\n")
        for r in rows:
            vals = [r.get(k, "-") for k in fields]
            f.write("| " + " | ".join(str(v).replace("|", "\\|") for v in vals) + " |\n")
    passed = sum(r["status"] == "PASS" for r in rows)
    validation_passed = sum(r["instance"] == "H" and r["status"] == "PASS" for r in rows)
    algorithm_passed = passed - validation_passed
    print(f"\nSummary: {algorithm_passed} algorithm runs passed; "
          f"{validation_passed} H validation passed; {failures} failed.")
    print(f"Saved: {csv_path.relative_to(HERE)}")
    print(f"Saved: {(RESULTS / 'summary.md').relative_to(HERE)}")
    return 1 if failures else 0

if __name__ == "__main__":
    raise SystemExit(main())
