# SARVault

[![CI](https://github.com/Knusperftw/SARVault/actions/workflows/ci.yml/badge.svg)](https://github.com/Knusperftw/SARVault/actions/workflows/ci.yml)

📚 **[Data docs & lineage](https://knusperftw.github.io/SARVault/)** (dbt docs, published from CI)

A reproducible, layered data warehouse over the public **ChEMBL** bioactivity
database, scoped to **cytotoxic / tubulin-targeting compounds** - the chemical
space behind antibody–drug-conjugate (ADC) payloads and classical
chemotherapeutics.

This is a data-engineering portfolio project. See [`SPEC.md`](./SPEC.md) for the
full design.

## Why this project

It demonstrates ingestion from a real external REST API (pagination,
idempotency, provenance), a medallion-style dbt transformation layer, a
documented dimensional model, and a warehouse that runs on **DuckDB**
(local / CI) and is **Snowflake-ready** - the same dbt models build against a
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

# 2. compute ECFP4 fingerprints + Murcko scaffolds into raw/ (RDKit)
python -m extract.cheminfo

# 3. project the fingerprints into a 2-D UMAP embedding (chemical space)
python -m extract.embedding

# 4. build the warehouse (staging -> marts -> analytics) + run dbt tests
dbt build --project-dir dbt --profiles-dir dbt/profiles

# 5. launch the dashboard
streamlit run dashboard/app.py
```

Or drive the whole lineage through Dagster (see [`docs/ORCHESTRATION.md`](./docs/ORCHESTRATION.md)):

```bash
dagster dev            # UI + lineage graph at http://localhost:3000
```

## Project status

Core pipeline **built end to end**: config-driven ChEMBL extract with provenance;
dbt medallion (staging → star-schema marts → analytical marts) on DuckDB; an RDKit
cheminfo stage deriving ECFP4 fingerprints and Bemis–Murcko scaffolds
(`dim_scaffold`, `mart_compound_fingerprint`) as the substrate for similarity,
substructure and scaffold analytics; UniChem and PDBe cross-reference enrichment
with an embedded 3D co-crystal viewer; and a
Streamlit dashboard (deployed on Streamlit Community Cloud). The full lineage is
orchestrated as a **Dagster asset graph** with dbt tests surfaced as asset checks
(see [`docs/ORCHESTRATION.md`](./docs/ORCHESTRATION.md)). The warehouse is
Snowflake-ready via a second dbt profile, not yet run against a live account. The
activity-cliff floor is set from the real-data profile (`scripts/profile_sar.py`)
via dbt vars (`cliff_min_tanimoto`, `cliff_min_delta_pchembl`).

## Exploring structure

The **Compound Library** page turns the ECFP4 fingerprints and Murcko scaffolds
into interactive structure search:

- **Structural analogs** - open any compound to see its nearest neighbours by ECFP4
  Tanimoto similarity, each annotated with its best potency and the potency delta to
  the query (a positive Δ flags a more potent close analog - a lead for the series).
- **Substructure filter** - enter a SMARTS query in the sidebar to keep only compounds
  containing that motif, with the matching atoms highlighted in the 2D depiction. For
  example, `c1ccccc1` keeps benzene-bearing compounds, `C(=O)N` keeps amides, and
  `[#7]` keeps anything with a nitrogen. An invalid pattern is flagged and ignored.
- **Activity cliffs** - the **Activity Cliffs** page surfaces pairs of structurally
  similar compounds (high ECFP4 Tanimoto) whose potency differs sharply on the same
  target - the sharpest signal in SAR. Pairs are ranked by SALI (= |Δ pChEMBL| /
  (1 − Tanimoto)) and shown side by side; similarity and Δ-potency thresholds are
  adjustable, and identical-2D-fingerprint pairs (stereo/tautomer/replicate) are
  flagged rather than mistaken for structural cliffs.
- **Chemical series** - the **Chemical Series** page groups compounds by their
  Bemis-Murcko scaffold, showing each series' size, potency spread and target reach,
  and drilling into the shared scaffold and its member compounds.
- **Structural embedding** - the **Chemical Space** page can project every compound
  into a 2-D UMAP of its ECFP4 fingerprint (proximity ≈ structural similarity),
  coloured by potency, approval class or scaffold series, so it's visible where
  potency concentrates in chemical space. Coordinates are precomputed once in the
  pipeline (fixed seed) and the physicochemical property scatter remains a toggle.

To characterise a real build (survival through the filters, scaffold series vs.
singletons, per-target pair budget, and the activity-cliff count across a
Tanimoto × Δ-pChEMBL grid), run the profiler against the warehouse:

```bash
python -m scripts.profile_sar --db warehouse.duckdb --out docs/DATA_PROFILE.md
```

## Change tracking across ChEMBL releases

The warehouse models *how the data changes over time*, not just its current state:

- **Incremental fact.** `fact_activity` is a dbt **incremental** model keyed on
  `activity_id`. ChEMBL activities are append-only with a stable id, so a refresh only
  scans activities not already loaded: re-running an unchanged release is a no-op, and a
  new release inserts just its delta (idempotent via `delete+insert` on the key).
- **SCD2 status history.** A dbt **snapshot** (`compound_status`) records a Type-2
  history of each compound's development status (`max_phase`, `is_approved_drug`,
  `pref_name`), keyed on the ChEMBL natural key `molecule_chembl_id`. When a payload
  advances — e.g. from a research compound to an approved drug — the snapshot closes the
  old validity window and opens a new one, so the warehouse can answer *when* a
  compound's approval status changed across releases.

The two behaviours are proven end to end against two pinned "releases": the v1 fixture
and a documented v2 delta (`tests/fixtures/build_raw_v2.py`, in which `CHEMBLM2` advances
`max_phase` 2 → 4 and two new activities appear). `tests/test_snapshot.py` builds both on
one DuckDB file and asserts the no-op re-run, the delta insert, and the second SCD2
window. To run the demo locally:

```bash
python tests/fixtures/build_raw_v2.py          # (re)generate the v2 delta fixture
pytest tests/test_snapshot.py -q               # incremental + SCD2 across two releases
```

## Data observability

Three layers of "is the data right *today*?" run in CI:

- **Source freshness.** Every raw extract stamps an ISO-8601 `_fetch_ts`; `dbt source
  freshness` measures the staleness of all nine raw sources against a batch-refresh
  cadence, so a forgotten refresh surfaces as a warning rather than silently serving old
  data.
- **Anomaly guards.** Three singular tests (`dbt/tests/assert_*.sql`) fail the build on a
  broken pipeline rather than bad numbers: a per-layer **row-count floor** (catches a
  silently-empty extract), a **null-rate ceiling** on `canonical_smiles` (catches a
  degraded structure layer), and a **pChEMBL mean band** (catches a unit error or bad
  merge). Thresholds are dbt vars, sized for CI and raised in production.
  `tests/test_observability.py` proves they fire by building against a seeded bad fixture
  (`tests/fixtures/build_raw_bad.py`) and asserting each one fails.
- **Published docs + lineage.** CI generates the dbt docs site (model descriptions and
  the full lineage graph) and deploys it to GitHub Pages:
  **[knusperftw.github.io/SARVault](https://knusperftw.github.io/SARVault/)**.

## Data provenance & license

Bioactivity data originates from **ChEMBL** (EMBL-EBI), pinned to a specific
release via [`config/target_set.yml`](./config/target_set.yml) (`chembl_version`).
ChEMBL data is released under a Creative Commons Attribution-ShareAlike license —
see [`LICENSE-DATA.md`](./LICENSE-DATA.md). No bulk ChEMBL data is committed to
this repository. Project code is licensed under MIT (see [`LICENSE`](./LICENSE)).
