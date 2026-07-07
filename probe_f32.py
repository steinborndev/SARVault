"""F3.2 ground-truth probe (stdlib only).

Run from the SARVault repo root:
    python probe_f32.py

Optionally point it at your built warehouse to also check which reference
payloads are already present:
    DUCKDB_PATH=warehouse.duckdb python probe_f32.py

Paste the whole output back. It resolves reference-payload ChEMBL ids, checks
warehouse presence, and reports the cellular-assay data shape + volume so the
F3.2 extract can be built against the real data rather than guessed.
"""

import json
import os
import urllib.parse
import urllib.request

BASE = "https://www.ebi.ac.uk/chembl/api/data"

# Reference ADC payloads to anchor (name -> intended payload_class).
REFERENCE_PAYLOADS = [
    ("CAMPTOTHECIN", "topo1_inhibitor"),
    ("SN-38", "topo1_inhibitor"),
    ("EXATECAN", "topo1_inhibitor"),
    ("TOPOTECAN", "topo1_inhibitor"),
    ("IRINOTECAN", "topo1_inhibitor"),
    ("BELOTECAN", "topo1_inhibitor"),
    ("MONOMETHYL AURISTATIN E", "tubulin_inhibitor"),
    ("MONOMETHYL AURISTATIN F", "tubulin_inhibitor"),
    ("MERTANSINE", "tubulin_inhibitor"),
    ("MAYTANSINE", "tubulin_inhibitor"),
]


def _get(path):
    url = f"{BASE}/{path}"
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=60) as resp:
        return json.load(resp)


def _count(path):
    return _get(path).get("page_meta", {}).get("total_count", 0)


def resolve(name):
    """Resolve a name to a molecule_chembl_id via exact pref_name, then synonym search."""
    q = urllib.parse.quote(name)
    exact = _get(f"molecule?pref_name__iexact={q}&format=json&limit=5").get("molecules", [])
    if exact:
        return [(m["molecule_chembl_id"], m.get("pref_name")) for m in exact]
    syn = _get(f"molecule/search?q={q}&format=json&limit=3").get("molecules", [])
    return [(m["molecule_chembl_id"], m.get("pref_name")) for m in syn]


def cellular_shape(chembl_id):
    """Counts of cellular-ish endpoints + one sample record's fields for a molecule."""
    gi50 = _count(f"activity?molecule_chembl_id={chembl_id}&standard_type=GI50&limit=1")
    ic50 = _count(f"activity?molecule_chembl_id={chembl_id}&standard_type=IC50&limit=1")
    ic50_f = _count(
        f"activity?molecule_chembl_id={chembl_id}&standard_type=IC50&assay_type=F&limit=1"
    )
    sample = _get(
        f"activity?molecule_chembl_id={chembl_id}&standard_type=GI50&limit=1"
    ).get("activities", [])
    keys = sorted(sample[0].keys()) if sample else []
    slim = (
        {
            k: sample[0].get(k)
            for k in (
                "assay_type", "target_type", "target_chembl_id", "target_pref_name",
                "target_organism", "assay_chembl_id", "standard_type",
                "standard_value", "standard_units", "pchembl_value",
            )
        }
        if sample
        else {}
    )
    return {"GI50": gi50, "IC50": ic50, "IC50_functional": ic50_f, "sample": slim, "all_keys": keys}


def warehouse_presence(ids):
    db = os.environ.get("DUCKDB_PATH")
    if not db or not os.path.exists(db):
        return None
    try:
        import duckdb
    except ImportError:
        return "duckdb not importable in this env"
    con = duckdb.connect(db, read_only=True)
    idset = ",".join(f"'{i}'" for i in ids if i)
    rows = con.execute(
        f"select molecule_chembl_id, pref_name from main_marts.dim_compound "
        f"where molecule_chembl_id in ({idset})"
    ).fetchall()
    return rows


def main():
    print("=== reference-payload id resolution ===")
    resolved = {}
    for name, cls in REFERENCE_PAYLOADS:
        try:
            hits = resolve(name)
        except Exception as exc:
            print(f"{name:26} ERROR {exc}")
            continue
        top = hits[0] if hits else (None, None)
        resolved[name] = top[0]
        alt = "" if len(hits) <= 1 else f"  (alts: {[h[0] for h in hits[1:]]})"
        print(f"{name:26} -> {top[0]}  [{top[1]}]  class={cls}{alt}")

    print("\n=== warehouse presence (dim_compound) ===")
    pres = warehouse_presence(list(resolved.values()))
    if pres is None:
        print("set DUCKDB_PATH to your warehouse to check presence")
    else:
        print(pres if pres else "none of the resolved ids are in the warehouse")

    print("\n=== cellular-assay shape + volume (first resolvable id per class) ===")
    seen_classes = set()
    for name, cls in REFERENCE_PAYLOADS:
        cid = resolved.get(name)
        if not cid or cls in seen_classes:
            continue
        seen_classes.add(cls)
        try:
            shape = cellular_shape(cid)
        except Exception as exc:
            print(f"{name} ({cid}) ERROR {exc}")
            continue
        print(f"\n{name} ({cid}, {cls}):")
        print(f"  GI50={shape['GI50']}  IC50={shape['IC50']}  IC50@assay_type=F={shape['IC50_functional']}")
        print(f"  sample record: {json.dumps(shape['sample'])}")
        print(f"  all record keys: {shape['all_keys']}")


if __name__ == "__main__":
    main()
