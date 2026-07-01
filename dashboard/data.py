"""Read-only data access for the SARVault dashboard (queries the DuckDB warehouse)."""

import os

import duckdb
import pandas as pd

DEFAULT_DB = os.environ.get("DUCKDB_PATH", "warehouse.duckdb")


def connect(db_path: str | None = None) -> duckdb.DuckDBPyConnection:
    """Open a read-only connection to the warehouse."""
    return duckdb.connect(str(db_path or DEFAULT_DB), read_only=True)


def load_target_sar(con: duckdb.DuckDBPyConnection) -> pd.DataFrame:
    return con.execute("select * from main_analytics.mart_target_sar").df()


def load_selectivity(con: duckdb.DuckDBPyConnection) -> pd.DataFrame:
    return con.execute("select * from main_analytics.mart_compound_selectivity").df()


def load_chemical_space(con: duckdb.DuckDBPyConnection) -> pd.DataFrame:
    return con.execute("select * from main_analytics.mart_chemical_space").df()


def load_compound_catalog(con: duckdb.DuckDBPyConnection) -> pd.DataFrame:
    return con.execute("select * from main_analytics.mart_compound_catalog").df()


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
    """Co-crystal PDB entries for one compound (its resolved PDBe structures)."""
    return con.execute(
        """
        select distinct ligand_code, pdb_id
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
    rows = [
        {"table": name, "rows": con.execute(f"select count(*) from {ref}").fetchone()[0]}
        for name, ref in tables.items()
    ]
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
