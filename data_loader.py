"""Load and NORMALISE the DLP instances from assignment_data.py, whatever its layout.

run_all.py, the tests and the experiments all go through load_instances(), so a different
assignment_data.py (for example the professor's) works without editing any solver.

Where the file is looked for (first hit wins; a notice is printed if several different ones exist):
    1. --data PATH or the environment variable ASSIGNMENT_DATA
    2. the current directory     3. the folder of this file     4. its parent folder
Layouts understood (labels such as A, B, C, C2, D ... are case-insensitive; 'Instance A', 'C_2' work):
    INSTANCES = {"A": {"p":..,"g":..,"h":..,"order":..}, ...}      (any dict-of-dicts, any variable name)
    INSTANCES = {"A": (p, g, h)  or  (p, g, h, order), ...}
    one module variable per instance:  A = {...}   INSTANCE_A = (...)   dlp_c2 = {...}
    one module variable per field:     A_p, A_g, A_h, A_order ...
    a list of dicts, each with a 'label' / 'name' / 'id' key
    a function get_instances() / load_instances() / instances() returning any of the above
Field names accepted (case-insensitive): p|prime|modulus, g|generator|base|alpha, h|target|beta|y,
order|n|group_order|ord, factors|factorization|factorisation.
Only p, g, h are required: the order defaults to p-1 and is replaced by the TRUE order of g, and
the factorisation is computed (Pollard rho) when it is not supplied. Secrets are never read here.
"""
import hashlib
import importlib.util
import os
import re
import sys

from common import factorize, element_order, is_probable_prime


class DataFormatError(RuntimeError):
    pass


_P = {"p", "prime", "modulus", "mod", "pr"}
_G = {"g", "generator", "gen", "base", "alpha"}
_H = {"h", "target", "beta", "y"}
_N = {"order", "n", "grouporder", "ord", "orderg", "ordg", "subgrouporder", "ordern"}
_F = {"factors", "factorization", "factorisation", "fac", "orderfactors", "primefactors", "nfactors"}
_HINTS = {"B_ic": {"bic", "bindexcalculus"}, "B_ls": {"bls"}, "M_ls": {"mls"}}
_SKIP_ATTR = {"SHA_VECTORS", "MESSAGES", "TEST_SECRETS", "SECRETS"}


def _key(k):
    return re.sub(r"[^a-z0-9]", "", str(k).lower())


def canon_label(s):
    """'Instance A' / 'inst_c2' / 'DLP-C_2' -> 'A' / 'C2' / 'C2'; None if it is not a label."""
    t = re.sub(r"[^a-z0-9]", "", str(s).lower())
    t = re.sub(r"^(instance|inst|dlp|problem|case|test)", "", t)
    m = re.fullmatch(r"([a-z])(\d{0,2})", t)
    return (m.group(1).upper() + m.group(2)) if m else None


def _to_int(v):
    if isinstance(v, bool):
        raise ValueError
    if isinstance(v, int):
        return v
    if isinstance(v, str):
        return int(v.strip().replace("_", "").replace(",", ""), 0)
    raise ValueError


def _mapping_of(obj):
    if isinstance(obj, dict):
        return obj
    if hasattr(obj, "_asdict"):
        return obj._asdict()
    if hasattr(obj, "__dict__") and not isinstance(obj, type):
        return vars(obj)
    try:                                         # MappingProxyType and friends
        return dict(obj.items())
    except Exception:
        return None


def _normalise_value(v):
    """dict / namedtuple / object / sequence -> {'p','g','h',...} or None."""
    m = _mapping_of(v)
    if m is not None:
        f = {}
        for k, val in m.items():
            kk = _key(k)
            for name, aliases in (("p", _P), ("g", _G), ("h", _H), ("order", _N)):
                if kk in aliases and name not in f:
                    try:
                        f[name] = _to_int(val)
                    except (ValueError, TypeError):
                        pass
            if kk in _F and "factors" not in f and isinstance(val, dict):
                try:
                    f["factors"] = {_to_int(a): _to_int(b) for a, b in val.items()}
                except (ValueError, TypeError):
                    pass
            for hint, aliases in _HINTS.items():
                if kk in aliases:
                    try:
                        f[hint] = _to_int(val)
                    except (ValueError, TypeError):
                        pass
        if {"p", "g", "h"} <= set(f):
            return f
        for val in m.values():                    # one level of nesting: {"params": {...}}
            if isinstance(val, (dict,)) or hasattr(val, "_asdict"):
                sub = _normalise_value(val)
                if sub:
                    return sub
        return None
    if isinstance(v, (tuple, list)) and 3 <= len(v) <= 6:
        try:
            nums = [_to_int(x) for x in v[:4] if not isinstance(x, (dict, list, tuple))]
        except (ValueError, TypeError):
            return None
        if len(nums) >= 3:
            first3 = nums[:3]
            p = max(first3)                       # the modulus is the largest of (p, g, h)
            g, h = [x for i, x in enumerate(first3) if i != first3.index(p)]
            f = {"p": p, "g": g, "h": h}
            if len(nums) == 4 and 1 < nums[3] <= p:
                f["order"] = nums[3]
            return f
    return None


def _scan(mod):
    """Return {label: raw_fields} found in the module, trying each supported layout."""
    for hook in ("get_instances", "load_instances", "instances"):
        fn = getattr(mod, hook, None)
        if callable(fn):
            try:
                found = _from_container(fn())
                if found:
                    return found
            except Exception:
                pass
    best = {}
    names = [n for n in dir(mod) if not n.startswith("_") and n not in _SKIP_ATTR]
    for n in names:                               # (a) containers: dict / list of instances
        found = _from_container(getattr(mod, n))
        if len(found) > len(best):
            best = found
    if best:
        return best
    found = {}                                    # (b) one variable per instance
    for n in names:
        lab = canon_label(n)
        if lab:
            nv = _normalise_value(getattr(mod, n))
            if nv:
                found[lab] = nv
    if found:
        return found
    fields = {}                                   # (c) one variable per field: A_p, g_C2 ...
    for n in names:
        m = re.fullmatch(r"(.+?)[_]+([A-Za-z0-9]+)", n)
        if not m:
            continue
        for lab_part, fld_part in ((m.group(1), m.group(2)), (m.group(2), m.group(1))):
            lab = canon_label(lab_part)
            fk = _key(fld_part)
            if lab and any(fk in al for al in (_P, _G, _H, _N)):
                try:
                    fields.setdefault(lab, {})[fk] = _to_int(getattr(mod, n))
                except (ValueError, TypeError):
                    pass
    out = {}
    for lab, d in fields.items():
        nv = _normalise_value(d)
        if nv:
            out[lab] = nv
    return out


def _from_container(obj):
    """{label: fields} from a container of instances: a mapping keyed by label, or a list of
    mappings that each carry their own label ('label' / 'name' / 'id' / 'instance' / 'case')."""
    found = {}
    if isinstance(obj, dict) or type(obj).__name__ == "mappingproxy":
        for key, value in obj.items():
            lab = canon_label(key)
            fields = _normalise_value(value) if lab else None
            if fields:
                found[lab] = fields
    elif isinstance(obj, (list, tuple)) and obj:
        for item in obj:
            mapping = _mapping_of(item)
            if mapping is None:
                return {}                         # not a list of instances (e.g. a tuple of numbers)
            for k in ("label", "name", "id", "instance", "case"):
                if k in mapping and canon_label(mapping[k]):
                    fields = _normalise_value(item)
                    if fields:
                        found[canon_label(mapping[k])] = fields
                    break
    return found


# --------------------------------------------------------------------------- file handling
def find_data_file(explicit=None, quiet=False):
    here = os.path.dirname(os.path.abspath(__file__))
    explicit = explicit or os.environ.get("ASSIGNMENT_DATA")
    if explicit:
        if not os.path.isfile(explicit):
            raise DataFormatError(f"data file not found: {explicit}")
        return os.path.abspath(explicit)
    cands = []
    for d in (os.getcwd(), here, os.path.dirname(here)):
        f = os.path.join(d, "assignment_data.py")
        if os.path.isfile(f) and os.path.abspath(f) not in cands:
            cands.append(os.path.abspath(f))
    if not cands:
        raise DataFormatError("assignment_data.py not found (looked in the current folder, "
                              f"{here} and its parent). Use --data PATH.")
    if len(cands) > 1 and not quiet:
        digests = {hashlib.md5(open(c, "rb").read()).hexdigest() for c in cands}
        if len(digests) > 1:
            print("NOTICE: several different assignment_data.py files exist; using the first:\n  "
                  + "\n  ".join(cands), file=sys.stderr)
    return cands[0]


def load_module(path):
    spec = importlib.util.spec_from_file_location("assignment_data", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["assignment_data"] = mod
    spec.loader.exec_module(mod)
    return mod


def _finish(label, f):
    """Add the derived fields every solver needs and sanity-check the instance."""
    p, g, h = f["p"], f["g"], f["h"] % f["p"]
    inst = {"label": label, "p": p, "g": g % p, "h": h, "notes": []}
    for k in ("B_ic", "B_ls", "M_ls"):
        if k in f:
            inst[k] = f[k]
    inst["p_is_prime"] = is_probable_prime(p)
    if not inst["p_is_prime"]:
        inst["notes"].append("p is not prime")
    given = f.get("order") or (p - 1)
    if pow(inst["g"], given, p) != 1:
        inst["notes"].append(f"supplied order {given} is not a multiple of ord(g); using p-1")
        given = p - 1
    inst["given_order"] = given
    facs = f.get("factors")
    prod = 1
    if facs:
        for q, e in facs.items():
            prod *= q ** e
    try:
        fac_given = facs if (facs and prod == given and all(is_probable_prime(q) for q in facs)) \
            else factorize(given)
        order = element_order(inst["g"], p, given, fac_given)
        inst["order_exact"] = True
        inst["factors"] = factorize(order) if order != given else dict(fac_given)
    except ValueError:                           # could not factor: keep the supplied order
        order, inst["order_exact"] = given, False
        inst["factors"] = None
        inst["notes"].append("order could not be factored in time")
    inst["order"] = order
    inst["generator"] = order == p - 1
    inst["safe_prime"] = bool(inst["p_is_prime"] and is_probable_prime((p - 1) // 2) and (p - 1) % 2 == 0)
    inst["h_in_subgroup"] = pow(inst["h"], order, p) == 1
    if not inst["h_in_subgroup"]:
        inst["notes"].append("h is not in the subgroup generated by g: no solution exists")
    return inst


def load_instances(path=None, quiet=False):
    """Return ({label: instance}, path_used). Raises DataFormatError with guidance if nothing fits."""
    path = find_data_file(path, quiet)
    try:
        mod = load_module(path)
    except Exception as exc:                     # syntax error, missing import, ... in the data file itself
        raise DataFormatError(f"could not import {path}: {type(exc).__name__}: {exc}") from exc
    try:
        raw = _scan(mod)
    except Exception as exc:
        raise DataFormatError(f"could not interpret {path}: {type(exc).__name__}: {exc}") from exc
    if not raw:
        shown = ", ".join(n for n in dir(mod) if not n.startswith("_"))[:300]
        raise DataFormatError(
            f"could not find DLP instances in {path}.\nPublic names in the file: {shown}\n"
            "Expected a dict such as INSTANCES = {'A': {'p':..., 'g':..., 'h':...}, ...} "
            "(see the docstring of data_loader.py for every supported layout).")
    ordered = dict(sorted(raw.items(), key=lambda kv: (kv[0][0], len(kv[0]), kv[0])))
    return {lab: _finish(lab, f) for lab, f in ordered.items()}, path
