"""Profile the built SARVault warehouse to characterise the SAR / cheminfo layer.

The point of this script is to replace guesses with numbers before building the
activity-cliff layer (F1.4). Run it against a warehouse built from the *real*
ChEMBL extract (with the cheminfo stage), not the CI fixtures, and it reports:

  1. Volume & survival        — rows per layer, cheminfo coverage.
  2. Scaffold distribution    — singletons vs. real chemical series.
  3. Per-target pair budget   — the pairwise self-join cost F1.4 will pay.
  4. Activity-cliff matrix    — pair counts across a Tanimoto × Δ-pChEMBL grid,
                                so the thresholds are chosen from the data.
  5. Top cliffs by SALI       — concrete examples to sanity-check the concept.

Tanimoto is computed with integer popcount over the ``ecfp4_hex`` fingerprint
(the same definition the dashboard uses), reconstructing each fingerprint once.

Prerequisite: the warehouse must contain ``mart_compound_fingerprint`` and
``dim_scaffold`` (added in the cheminfo PR). A warehouse built before that layer
existed will not have them — rebuild locally first:

    python -m extract.run                 # real ChEMBL slice (or reuse raw/)
    python -m extract.cheminfo            # ECFP4 + Murcko scaffolds
    DUCKDB_PATH=warehouse.duckdb dbt build --project-dir dbt --profiles-dir dbt/profiles

Usage:
    python -m scripts.profile_sar --db warehouse.duckdb --out docs/DATA_PROFILE.md
"""

import argparse
import heapq
import sys
from itertools import combinations
from pathlib import Path

import duckdb

DEFAULT_TAN_THRESHOLDS = (0.70, 0.75, 0.80, 0.85, 0.90)
DEFAULT_DELTA_THRESHOLDS = (1.0, 1.5, 2.0, 3.0)


# --- fingerprint similarity (integer popcount over the hex fingerprint) ---------
def tanimoto(a_int: int, pa: int, b_int: int, pb: int) -> float:
    """Tanimoto from precomputed fingerprint ints and their popcounts."""
    inter = (a_int & b_int).bit_count()
    union = pa + pb - inter
    return inter / union if union else 0.0


# --- data access ----------------------------------------------------------------
def _has_relation(con, relation: str) -> bool:
    try:
        con.execute(f"select 1 from {relation} limit 1")
        return True
    except Exception:
        return False


def require_cheminfo(con) -> None:
    missing = [
        r
        for r in ("main_analytics.mart_compound_fingerprint", "main_marts.dim_scaffold")
        if not _has_relation(con, r)
    ]
    if missing:
        sys.exit(
            "This warehouse has no cheminfo layer ("
            + ", ".join(missing)
            + ").\nRebuild with `python -m extract.cheminfo` + `dbt build` first "
            "(see the module docstring)."
        )


def load_volume(con) -> dict:
    q = {
        "compounds (dim_compound)": "select count(*) from main_marts.dim_compound",
        "targets (dim_target)": "select count(*) from main_marts.dim_target",
        "assays (dim_assay)": "select count(*) from main_marts.dim_assay",
        "activities (fact_activity)": "select count(*) from main_marts.fact_activity",
        "compound-target pairs (mart_target_sar)": "select count(*) from main_analytics.mart_target_sar",
        "fingerprinted compounds": "select count(*) from main_analytics.mart_compound_fingerprint",
        "distinct scaffolds (dim_scaffold)": "select count(*) from main_marts.dim_scaffold",
    }
    return {label: con.execute(sql).fetchone()[0] for label, sql in q.items()}


def load_scaffold_sizes(con):
    return con.execute(
        "select n_compounds from main_marts.dim_scaffold order by n_compounds desc"
    ).df()["n_compounds"].tolist()


def load_target_compounds(con):
    """Per target: rows of (compound_key, molecule_chembl_id, ecfp4_hex, median_pchembl).

    Joins the per-target potency summary to the fingerprint mart; compounds without
    a fingerprint (unparseable SMILES) are dropped and counted separately.
    """
    rows = con.execute(
        """
        select
            s.target_key,
            t.pref_name              as target,
            s.compound_key,
            s.molecule_chembl_id,
            f.ecfp4_hex,
            s.median_pchembl
        from main_analytics.mart_target_sar s
        join main_marts.dim_target t on s.target_key = t.target_key
        left join main_analytics.mart_compound_fingerprint f
               on s.compound_key = f.compound_key
        """
    ).df()
    n_missing_fp = int(rows["ecfp4_hex"].isna().sum())
    rows = rows.dropna(subset=["ecfp4_hex", "median_pchembl"])

    per_target = {}
    for rec in rows.itertuples(index=False):
        per_target.setdefault((rec.target_key, rec.target), []).append(
            (rec.compound_key, rec.molecule_chembl_id, rec.ecfp4_hex, float(rec.median_pchembl))
        )
    return per_target, n_missing_fp


# --- cliff profiling ------------------------------------------------------------
def profile_cliffs(per_target, tan_thresholds, delta_thresholds, top_k, max_per_target):
    """Walk per-target compound pairs, tallying the cliff grid and top-SALI cliffs.

    Returns (grid, top_cliffs, pair_budget) where grid[(t, d)] is the pair count
    with Tanimoto >= t and |Δ pChEMBL| >= d, top_cliffs is a SALI-ranked list, and
    pair_budget lists (target, n_compounds, n_pairs) for the self-join cost.
    """
    grid = {(t, d): 0 for t in tan_thresholds for d in delta_thresholds}
    min_tan = min(tan_thresholds)
    min_delta = min(delta_thresholds)
    top = []  # min-heap of (sali, target, id_a, id_b, tan, delta)
    pair_budget = []

    for (_tkey, target), members in sorted(per_target.items(), key=lambda kv: -len(kv[1])):
        if max_per_target and len(members) > max_per_target:
            # Keep the most potent compounds (most SAR-relevant) when capping.
            members = sorted(members, key=lambda m: -m[3])[:max_per_target]
        n = len(members)
        pair_budget.append((target, n, n * (n - 1) // 2))
        if n < 2:
            continue

        # Precompute fingerprint ints + popcounts once per compound.
        prepped = [
            (ck, cid, (fi := int(hexv, 16)), fi.bit_count(), pot)
            for ck, cid, hexv, pot in members
        ]

        for (ka, ida, ia, pa, pota), (kb, idb, ib, pb, potb) in combinations(prepped, 2):
            t = tanimoto(ia, pa, ib, pb)
            if t < min_tan:
                continue
            delta = abs(pota - potb)
            if delta < min_delta:
                continue
            for tt in tan_thresholds:
                if t < tt:
                    continue
                for dd in delta_thresholds:
                    if delta >= dd:
                        grid[(tt, dd)] += 1
            sali = delta / (1 - t) if t < 1.0 else float("inf")
            item = (sali, target, ida, idb, round(t, 3), round(delta, 2))
            if len(top) < top_k:
                heapq.heappush(top, item)
            elif sali > top[0][0]:
                heapq.heapreplace(top, item)

    top_cliffs = sorted(top, key=lambda x: -x[0])
    return grid, top_cliffs, pair_budget


# --- report rendering -----------------------------------------------------------
def _pct(part, whole):
    return f"{100 * part / whole:.1f}%" if whole else "—"


def render_report(volume, scaffold_sizes, n_missing_fp, grid, top_cliffs, pair_budget,
                  tan_thresholds, delta_thresholds) -> str:
    out = ["# SARVault — SAR data profile", ""]

    out += ["## 1. Volume & survival", ""]
    for label, n in volume.items():
        out.append(f"- {label}: **{n:,}**")
    fp = volume["fingerprinted compounds"]
    comp = volume["compounds (dim_compound)"]
    out.append(f"- fingerprint coverage: **{_pct(fp, comp)}** of compounds")
    out.append(f"- target-compound rows without a fingerprint (dropped): {n_missing_fp:,}")
    out.append("")

    out += ["## 2. Scaffold distribution", ""]
    n_scaf = len(scaffold_sizes)
    singletons = sum(1 for s in scaffold_sizes if s == 1)
    series = [s for s in scaffold_sizes if s >= 2]
    covered = sum(scaffold_sizes)
    out.append(f"- distinct scaffolds: **{n_scaf:,}**")
    out.append(f"- singletons (1 compound): {singletons:,} ({_pct(singletons, n_scaf)})")
    out.append(f"- real series (>=2 compounds): **{len(series):,}** ({_pct(len(series), n_scaf)})")
    if series:
        out.append(f"- largest series sizes: {sorted(series, reverse=True)[:10]}")
        out.append(f"- compounds sitting in a multi-member series: {sum(series):,} of {covered:,}")
    out.append("")

    out += ["## 3. Per-target pair budget (self-join cost)", ""]
    total_pairs = sum(p for _, _, p in pair_budget)
    out.append(f"- total pairwise comparisons across targets: **{total_pairs:,}**")
    out.append("")
    out.append("| target | compounds | pairs |")
    out.append("|---|---:|---:|")
    for target, n, p in pair_budget:
        out.append(f"| {target} | {n:,} | {p:,} |")
    out.append("")

    out += ["## 4. Activity-cliff matrix (pair counts)", ""]
    out.append("Pairs with Tanimoto >= row **and** |Δ pChEMBL| >= column.")
    out.append("")
    header = "| Tanimoto \\\\ Δ pChEMBL | " + " | ".join(f"≥{d}" for d in delta_thresholds) + " |"
    out.append(header)
    out.append("|---|" + "---:|" * len(delta_thresholds))
    for t in tan_thresholds:
        cells = " | ".join(f"{grid[(t, d)]:,}" for d in delta_thresholds)
        out.append(f"| ≥{t:.2f} | {cells} |")
    out.append("")

    out += ["## 5. Top cliffs by SALI", ""]
    if not top_cliffs:
        out.append("_No cliffs found at the lowest grid thresholds — see the matrix above._")
    else:
        out.append("SALI = |Δ pChEMBL| / (1 − Tanimoto). Higher = sharper cliff.")
        out.append("")
        out.append("| target | compound A | compound B | Tanimoto | Δ pChEMBL | SALI |")
        out.append("|---|---|---|---:|---:|---:|")
        for sali, target, ida, idb, t, d in top_cliffs:
            sali_s = "∞" if sali == float("inf") else f"{sali:.1f}"
            out.append(f"| {target} | {ida} | {idb} | {t} | {d} | {sali_s} |")
    out.append("")

    out += ["---", "",
            "_Generated by `scripts/profile_sar.py`. These numbers should drive the "
            "F1.4 activity-cliff thresholds and confirm the self-join stays tractable._"]
    return "\n".join(out) + "\n"


def main(argv=None):
    ap = argparse.ArgumentParser(description="Profile the SARVault warehouse for the SAR layer.")
    ap.add_argument("--db", default="warehouse.duckdb", help="path to the DuckDB warehouse")
    ap.add_argument("--out", default=None, help="write the markdown report to this path")
    ap.add_argument("--top-cliffs", type=int, default=25, help="how many top-SALI cliffs to list")
    ap.add_argument(
        "--max-per-target",
        type=int,
        default=0,
        help="cap compounds per target (0 = no cap); keeps the most potent when capping",
    )
    ap.add_argument("--tan-thresholds", default=",".join(str(t) for t in DEFAULT_TAN_THRESHOLDS))
    ap.add_argument("--delta-thresholds", default=",".join(str(d) for d in DEFAULT_DELTA_THRESHOLDS))
    args = ap.parse_args(argv)

    tan_thresholds = tuple(float(x) for x in args.tan_thresholds.split(","))
    delta_thresholds = tuple(float(x) for x in args.delta_thresholds.split(","))

    if not Path(args.db).exists():
        sys.exit(f"Warehouse not found: {args.db}")

    con = duckdb.connect(args.db, read_only=True)
    try:
        require_cheminfo(con)
        volume = load_volume(con)
        scaffold_sizes = load_scaffold_sizes(con)
        per_target, n_missing_fp = load_target_compounds(con)
    finally:
        con.close()

    grid, top_cliffs, pair_budget = profile_cliffs(
        per_target, tan_thresholds, delta_thresholds, args.top_cliffs, args.max_per_target
    )
    report = render_report(
        volume, scaffold_sizes, n_missing_fp, grid, top_cliffs, pair_budget,
        tan_thresholds, delta_thresholds,
    )

    print(report)
    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text(report)
        print(f"[written] {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
