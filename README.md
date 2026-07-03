# SARVault

[![CI](https://github.com/Knusperftw/SARVault/actions/workflows/ci.yml/badge.svg)](https://github.com/Knusperftw/SARVault/actions/workflows/ci.yml)

A reproducible, layered data warehouse over the public **ChEMBL** bioactivity
database, scoped to **cytotoxic / tubulin-targeting compounds** — the chemical
space behind antibody–drug-conjugate (ADC) payloads and classical
chemotherapeutics.

This is a data-engineering portfolio project. See [`SPEC.md`](./SPEC.md) for the
full design.

## Why this project

It demonstrates ingestion from a real external REST API (pagination,
idempotency, provenance), a medallion-style dbt transformation layer, a
documented dimensional model, and a warehouse that runs on **DuckDB**
(local / CI) and is **Snowflake-ready** — the same dbt models build against a
second `dbt-snowflake` profile (not yet run against a live Snowflake account).

It answers three questions about the payload chemical space: structure–activity
ranking per target, compound selectivity profiling, and chemical-space
characterization versus approved drugs.

## Architecture

```
ChEMBL REST API
   → raw (Parquet + provenance)
   → staging (dbt)
   → marts (star schema, dbt)
   → analytics (dbt)
   → Streamlit
                DuckDB  ⇄  Snowflake   (same models, swapped profile;
                                        Snowflake profile defined, not yet run)
```

## Quickstart

```bash
# install everything (add `orchestration` for the Dagster asset graph)
pip install -e ".[dev,extract,dbt,dashboard,orchestration]"

# lint + test
ruff check .
pytest -q
```

End-to-end, from the repo root:

```bash
# 1. extract the scoped ChEMBL slice into raw/
python -m extract.run

# 2. build the warehouse (staging -> marts -> analytics) + run dbt tests
dbt build --project-dir dbt --profiles-dir dbt/profiles

# 3. launch the dashboard
streamlit run dashboard/app.py
```

Or drive the whole lineage through Dagster (see [`docs/ORCHESTRATION.md`](./docs/ORCHESTRATION.md)):

```bash
dagster dev            # UI + lineage graph at http://localhost:3000
```

## Project status

Core pipeline **built end to end**: config-driven ChEMBL extract with provenance;
dbt medallion (staging → star-schema marts → analytical marts) on DuckDB; UniChem
and PDBe cross-reference enrichment with an embedded 3D co-crystal viewer; and a
Streamlit dashboard (deployed on Streamlit Community Cloud). The full lineage is
orchestrated as a **Dagster asset graph** with dbt tests surfaced as asset checks
(see [`docs/ORCHESTRATION.md`](./docs/ORCHESTRATION.md)). The warehouse is
Snowflake-ready via a second dbt profile, not yet run against a live account.

## Data provenance & license

Bioactivity data originates from **ChEMBL** (EMBL-EBI), pinned to a specific
release via [`config/target_set.yml`](./config/target_set.yml) (`chembl_version`).
ChEMBL data is released under a Creative Commons Attribution-ShareAlike license —
see [`LICENSE-DATA.md`](./LICENSE-DATA.md). No bulk ChEMBL data is committed to
this repository. Project code is licensed under MIT (see [`LICENSE`](./LICENSE)).
