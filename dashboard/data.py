"""Read-only data access for the SARVault dashboard (queries the DuckDB warehouse)."""

import os
from pathlib import Path

import duckdb
import pandas as pd

DEFAULT_DB = os.environ.get("DUCKDB_PATH", "warehouse.duckdb")


def connect(db_path: str | None = None) -> duckdb.DuckDBPyConnection:
    """Open a read-only connection to the warehouse."""
    return duckdb.connect(str(db_path or DEFAULT_DB), read_only=True)


def ensure_warehouse(
    db_path: str | None = None, url: str | None = None, token: str | None = None
) -> Path:
    """Ensure the warehouse file exists locally; fetch it from ``url`` if missing.

    On a fresh cloud deploy (e.g. Streamlit Community Cloud) the DuckDB file is
    not tracked in git, so it is downloaded once and cached in the container.
    Local runs that already have a warehouse are left untouched, and if no URL is
    provided the caller surfaces the usual "warehouse missing" error.

    ``url`` is a GitHub Release asset browser URL. If ``token`` is set, the asset
    is fetched authenticated via the GitHub API (required for private repos,
    whose browser download URLs 404 for anonymous requests); otherwise it is
    downloaded directly (public assets).
    """
    path = Path(db_path or DEFAULT_DB)
    if path.exists() or not url:
        return path
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".part")
    if token:
        _download_private_asset(url, token, tmp)
    else:
        import urllib.request

        urllib.request.urlretrieve(url, tmp)  # public asset, no auth
    tmp.replace(path)
    return path


def _parse_asset_url(browser_url: str) -> tuple[str, str, str, str]:
    """Split a release asset browser URL into (owner, repo, tag, asset_name)."""
    import re

    m = re.match(
        r"https://github\.com/([^/]+)/([^/]+)/releases/download/([^/]+)/(.+)$",
        browser_url,
    )
    if not m:
        raise ValueError(f"Unrecognized GitHub release asset URL: {browser_url}")
    return m.group(1), m.group(2), m.group(3), m.group(4)


def _download_private_asset(browser_url: str, token: str, dest: Path) -> None:
    """Download a private-repo release asset via the GitHub API, streaming to disk.

    Resolves the asset id for {tag}/{name}, then requests the API asset endpoint
    with ``Accept: application/octet-stream``. GitHub answers with a 302 to a
    signed storage URL on another host; ``requests`` drops the Authorization
    header on that cross-host redirect (which the storage backend requires).
    """
    import requests

    owner, repo, tag, name = _parse_asset_url(browser_url)
    api = f"https://api.github.com/repos/{owner}/{repo}/releases"
    auth = {"Authorization": f"Bearer {token}"}
    rel = requests.get(
        f"{api}/tags/{tag}",
        headers={**auth, "Accept": "application/vnd.github+json"},
        timeout=30,
    )
    rel.raise_for_status()
    asset = next((a for a in rel.json().get("assets", []) if a["name"] == name), None)
    if asset is None:
        raise FileNotFoundError(f"Asset {name!r} not found in release {tag!r} of {owner}/{repo}")
    with requests.get(
        asset["url"],
        headers={**auth, "Accept": "application/octet-stream"},
        stream=True,
        timeout=300,
    ) as resp:
        resp.raise_for_status()
        with open(dest, "wb") as fh:
            for chunk in resp.iter_content(chunk_size=1 << 20):
                fh.write(chunk)


def load_target_sar(con: duckdb.DuckDBPyConnection) -> pd.DataFrame:
    return con.execute("select * from main_analytics.mart_target_sar").df()


def load_selectivity(con: duckdb.DuckDBPyConnection) -> pd.DataFrame:
    return con.execute("select * from main_analytics.mart_compound_selectivity").df()


def load_chemical_space(con: duckdb.DuckDBPyConnection) -> pd.DataFrame:
    return con.execute("select * from main_analytics.mart_chemical_space").df()


def load_compound_catalog(con: duckdb.DuckDBPyConnection) -> pd.DataFrame:
    return con.execute("select * from main_analytics.mart_compound_catalog").df()


def load_fingerprints(con: duckdb.DuckDBPyConnection) -> pd.DataFrame:
    """Per-compound ECFP4 hex fingerprint + scaffold (empty if the mart is absent).

    A warehouse built before the cheminfo layer existed has no fingerprint mart;
    return an empty frame so the similarity/substructure features degrade quietly
    rather than breaking the whole page.
    """
    try:
        return con.execute(
            """
            select compound_key, molecule_chembl_id, ecfp4_hex,
                   murcko_scaffold_smiles, scaffold_key
            from main_analytics.mart_compound_fingerprint
            """
        ).df()
    except Exception:
        return pd.DataFrame(
            columns=[
                "compound_key",
                "molecule_chembl_id",
                "ecfp4_hex",
                "murcko_scaffold_smiles",
                "scaffold_key",
            ]
        )


def load_activity_cliffs(con: duckdb.DuckDBPyConnection) -> pd.DataFrame:
    """Activity-cliff pairs (empty if the mart is absent in an older warehouse)."""
    try:
        return con.execute(
            """
            select target_key, target_pref_name,
                   compound_key_a, compound_key_b,
                   molecule_chembl_id_a, molecule_chembl_id_b,
                   pchembl_a, pchembl_b, delta_pchembl, tanimoto,
                   is_identical_fp, same_scaffold, sali
            from main_analytics.mart_activity_cliff
            """
        ).df()
    except Exception:
        return pd.DataFrame(
            columns=[
                "target_key", "target_pref_name", "compound_key_a", "compound_key_b",
                "molecule_chembl_id_a", "molecule_chembl_id_b", "pchembl_a", "pchembl_b",
                "delta_pchembl", "tanimoto", "is_identical_fp", "same_scaffold", "sali",
            ]
        )


def load_chemical_series(con: duckdb.DuckDBPyConnection) -> pd.DataFrame:
    """Per-scaffold chemical series (empty if the mart is absent in an older warehouse)."""
    try:
        return con.execute("select * from main_analytics.mart_chemical_series").df()
    except Exception:
        return pd.DataFrame(
            columns=[
                "scaffold_key", "murcko_scaffold_smiles", "murcko_generic_smiles",
                "n_compounds", "n_measured_compounds", "n_targets", "median_pchembl",
                "max_pchembl", "min_pchembl", "pchembl_range", "top_compound",
            ]
        )


def scaffold_members(con: duckdb.DuckDBPyConnection, scaffold_key: int) -> pd.DataFrame:
    """Member compounds of one scaffold series, with structure and best potency."""
    return con.execute(
        """
        with potency as (
            select compound_key,
                   max(median_pchembl)        as best_pchembl,
                   count(distinct target_key) as n_targets
            from main_analytics.mart_target_sar
            group by compound_key
        )
        select
            f.molecule_chembl_id,
            c.pref_name,
            c.canonical_smiles,
            round(p.best_pchembl, 2) as best_pchembl,
            coalesce(p.n_targets, 0) as n_targets
        from main_analytics.mart_compound_fingerprint f
        join main_marts.dim_compound c on f.compound_key = c.compound_key
        left join potency p on f.compound_key = p.compound_key
        where f.scaffold_key = ?
        order by p.best_pchembl desc nulls last
        """,
        [scaffold_key],
    ).df()


def compound_row(con: duckdb.DuckDBPyConnection, molecule_chembl_id: str):
    """Full compound record for the shared detail card, by ChEMBL id.

    Prefers the catalog (has best_pchembl / selectivity_index / has_pdb); falls back
    to dim_compound for compounds absent from the catalog (e.g. a scaffold member
    with no measured activity), filling the catalog-only columns with nulls.
    """
    try:
        df = con.execute(
            "select * from main_analytics.mart_compound_catalog where molecule_chembl_id = ?",
            [molecule_chembl_id],
        ).df()
    except Exception:
        df = pd.DataFrame()  # catalog mart absent in this warehouse
    if len(df):
        return df.iloc[0]
    df = con.execute(
        "select *, cast(null as double) as selectivity_index "
        "from main_marts.dim_compound where molecule_chembl_id = ?",
        [molecule_chembl_id],
    ).df()
    return df.iloc[0] if len(df) else None


def compound_target_profile(con: duckdb.DuckDBPyConnection, compound_key: int) -> pd.DataFrame:
    """Per-target potency for one compound (its SAR fingerprint)."""
    return con.execute(
        """
        select
            t.pref_name               as target,
            round(s.median_pchembl, 2) as median_pchembl,
            round(s.max_pchembl, 2)    as max_pchembl,
            s.n_measurements,
            s.n_assays
        from main_analytics.mart_target_sar s
        join main_marts.dim_target t on s.target_key = t.target_key
        where s.compound_key = ?
        order by s.median_pchembl desc
        """,
        [compound_key],
    ).df()


def compound_xrefs(con: duckdb.DuckDBPyConnection, compound_key: int) -> pd.DataFrame:
    """One representative cross-reference per source, plus the total count per source."""
    return con.execute(
        """
        select
            display_name,
            min(xref_id)           as xref_id,
            arg_min(url, xref_id)  as url,
            count(*)               as n_refs
        from main_analytics.mart_compound_xref
        where compound_key = ?
        group by display_name
        order by display_name
        """,
        [compound_key],
    ).df()


def compound_pdb_entries(con: duckdb.DuckDBPyConnection, compound_key: int) -> pd.DataFrame:
    """Co-crystal PDB entries for one compound, enriched with entry metadata."""
    return con.execute(
        """
        select distinct ligand_code, pdb_id, title, method, year, resolution
        from main_analytics.mart_compound_pdb
        where compound_key = ?
        order by pdb_id
        """,
        [compound_key],
    ).df()


def list_target_names(con: duckdb.DuckDBPyConnection) -> list[str]:
    return (
        con.execute("select pref_name from main_marts.dim_target order by pref_name")
        .df()["pref_name"]
        .tolist()
    )


def target_summary(con: duckdb.DuckDBPyConnection) -> pd.DataFrame:
    """Per-target overview for the landing page."""
    return con.execute(
        """
        select
            t.pref_name                       as target,
            t.target_type,
            count(distinct s.compound_key)    as n_compounds,
            count(*)                          as n_pairs,
            round(median(s.median_pchembl),2) as median_pchembl,
            round(max(s.max_pchembl),2)       as best_pchembl
        from main_analytics.mart_target_sar s
        join main_marts.dim_target t on s.target_key = t.target_key
        group by t.pref_name, t.target_type
        order by n_pairs desc
        """
    ).df()


def confidence_distribution(con: duckdb.DuckDBPyConnection) -> pd.DataFrame:
    """Assay confidence distribution across the fact table."""
    return con.execute(
        """
        select
            da.confidence_score,
            count(*)                     as n_activities,
            count(distinct da.assay_key) as n_assays
        from main_marts.fact_activity f
        join main_marts.dim_assay da on f.assay_key = da.assay_key
        group by da.confidence_score
        order by da.confidence_score
        """
    ).df()


def standard_type_distribution(con: duckdb.DuckDBPyConnection) -> pd.DataFrame:
    """Activity and compound counts per ChEMBL standard (endpoint) type."""
    return con.execute(
        """
        select
            standard_type,
            count(*)                     as n_activities,
            count(distinct compound_key) as n_compounds
        from main_marts.fact_activity
        group by standard_type
        order by n_activities desc
        """
    ).df()


def layer_counts(con: duckdb.DuckDBPyConnection) -> pd.DataFrame:
    """Row counts per modelled layer, for the provenance view."""
    tables = {
        "stg_activities": "main_staging.stg_activities",
        "dim_compound": "main_marts.dim_compound",
        "dim_target": "main_marts.dim_target",
        "dim_assay": "main_marts.dim_assay",
        "fact_activity": "main_marts.fact_activity",
        "mart_target_sar": "main_analytics.mart_target_sar",
        "mart_compound_selectivity": "main_analytics.mart_compound_selectivity",
        "mart_chemical_space": "main_analytics.mart_chemical_space",
    }
    rows = []
    for name, ref in tables.items():
        try:
            n = con.execute(f"select count(*) from {ref}").fetchone()[0]
        except Exception:
            # A relation may be unreadable in a shipped snapshot (e.g. a staging
            # view over raw Parquet not present in the deployment); skip it rather
            # than failing the whole provenance page.
            continue
        rows.append({"table": name, "rows": n})
    return pd.DataFrame(rows)


def headline_metrics(con: duckdb.DuckDBPyConnection) -> dict:
    queries = {
        "compounds": "select count(*) from main_marts.dim_compound",
        "activities": "select count(*) from main_marts.fact_activity",
        "targets": "select count(*) from main_marts.dim_target",
        "approved": "select count(*) from main_marts.dim_compound where is_approved_drug",
    }
    return {key: con.execute(sql).fetchone()[0] for key, sql in queries.items()}


def pipeline_config() -> dict:
    """ChEMBL version / confidence floor / organism from the target-set config."""
    from extract.config import load_config

    cfg = load_config()
    return {
        "chembl_version": cfg.chembl_version,
        "min_confidence_score": cfg.activity.min_confidence_score,
        "organism": cfg.organism,
    }
