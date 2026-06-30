"""Dagster asset & job definitions for the end-to-end pipeline.

Materializes the medallion lineage (extract -> load_raw -> dbt build -> test).
Implemented in milestone M7 (feat/orchestration). Import-light by design.
"""

__all__: list[str] = []
