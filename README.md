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
documented dimensional model, and a warehouse that runs on both **DuckDB**
(local / CI) and **Snowflake** (cloud) from the same dbt models.

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
                DuckDB  ⇄  Snowflake   (same models, swapped profile)
```

## Quickstart

```bash
# install everything
pip install -e ".[dev,extract,dbt,dashboard]"

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

The pipeline stages are added across milestones M1–M7 — see [`SPEC.md`](./SPEC.md) §12.

## Project status

Milestone **M0 — scaffold**. Repository structure, container, CI, and the
target-set configuration are in place; ingestion and transformation logic land
in subsequent milestones.

## Data provenance & license

Bioactivity data originates from **ChEMBL** (EMBL-EBI), pinned to a specific
release via [`config/target_set.yml`](./config/target_set.yml) (`chembl_version`).
ChEMBL data is released under a Creative Commons Attribution-ShareAlike license —
see [`LICENSE-DATA.md`](./LICENSE-DATA.md). No bulk ChEMBL data is committed to
this repository. Project code is licensed under MIT (see [`LICENSE`](./LICENSE)).
