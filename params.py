"""Single source of truth for every tunable parameter used by run_all.py and the experiments.

The values are the STARTING PARAMETERS GIVEN IN THE ASSIGNMENT:
    index calculus : F -> B = 500,  G -> B = 2000
    linear sieve   : G -> B = 700,  M = 900
Other instances (for example the bonus instance H) get a heuristic B / the same B, M.
Hints stored in assignment_data.py (B_ic, B_ls, M_ls) are IGNORED unless --use-instance-params
is given, so run_all.py and the experiments can never silently disagree.
"""
IC_B = {"F": 500, "G": 2000}
LS_B, LS_M = 700, 900
BSGS_BYTES_PER_ENTRY = 120          # measured CPython dict cost per (int -> int) entry, rounded up


def ic_params(label, inst, use_instance_params=False):
    """Return (B, where_it_came_from)."""
    from index_calculus import suggest_B
    if use_instance_params and inst.get("B_ic"):
        return int(inst["B_ic"]), "data file"
    if label in IC_B:
        return IC_B[label], "assignment"
    return suggest_B(inst["p"]), "heuristic"


def ls_params(label, inst, use_instance_params=False):
    """Return ((B, M), where_it_came_from)."""
    if use_instance_params and inst.get("B_ls") and inst.get("M_ls"):
        return (int(inst["B_ls"]), int(inst["M_ls"])), "data file"
    return (LS_B, LS_M), "assignment"
