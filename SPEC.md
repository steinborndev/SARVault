# SPEC — ChEMBL Bioactivity SAR Warehouse

**Project:** `SARVault`
**Status:** Design draft (pre-implementation)
**Author:** Dr. Benjamin Steinborn
**Purpose:** Portfolio project demonstrating modern data-engineering competency on a real external life-science API, complementing the internal-R&D-focused `mAbVault` project.

---

## 1. Purpose & positioning

This project builds a reproducible, layered data warehouse over the public **ChEMBL** bioactivity database, scoped to **cytotoxic / tubulin-targeting compounds** — the chemical space underlying ADC payloads and classical chemotherapeutics.

It is designed to evidence:

- Ingestion from a real external REST API (pagination, idempotency, provenance) — not flat-file reads.
- A medallion-style transformation layer (raw → staging → marts) using dbt.
- A clean dimensional model with documented grain.
- A cloud-warehouse runtime target (Snowflake) alongside a local engine (DuckDB).
- Data-quality enforcement, automated testing, containerization, and CI.
- Domain-credible analytics (structure–activity relationships, selectivity profiling).

### Relationship to `mAbVault`

| Axis | `mAbVault` | This project |
|---|---|---|
| Data origin | Simulated internal R&D | Real external public API |
| Domain | mAb stability / CMC | Cytotoxic small-molecule SAR |
| Warehouse | DuckDB (local) | DuckDB **and** Snowflake (cloud) |
| Emphasis | Data *content* understanding | Modern DE *stack* breadth |

Together they signal: *"I understand pharma R&D data and I can build the modern stack against real sources."*

### Target roles

R&D Data Engineering (Roche) and CMC / analytical data roles (e.g. Tubulis). The payload-chemistry scope deliberately ties to the author's B.Sc. thesis (chondramides / microtubule targeting), Ph.D. (chemotherapeutic delivery), and ADC-payload relevance.

---

## 2. Analytical thesis

A pure "copy ChEMBL into a warehouse" pipeline has no narrative. This warehouse answers concrete questions about the payload chemical space:

1. **SAR ranking** — For a given target (e.g. tubulin / β-tubulin, topoisomerase, etc.), which compounds are most potent, and how is potency distributed across chemotypes?
2. **Selectivity profiling** — For compounds tested against multiple targets, what is the potency spread and selectivity index?
3. **Chemical-space characterization** — How are physicochemical properties (MW, logP, HBA/HBD, Ro5 violations) distributed across the payload set, and how do approved drugs sit within that space?

These three questions define the three analytical marts (Section 6.4).

---

## 3. Scope

### In scope

- A **config-driven target/mechanism set** (Section 6.0) defining the compound universe. Default: tubulin-targeting and selected cytotoxic-payload-relevant targets.
- Entities: molecules, targets, assays, activities (the ChEMBL core).
- Human-target focus by default (`organism = Homo sapiens`), configurable.
- Activity measurements with a standardized potency value (`pchembl_value`).
- Local (DuckDB) build as the primary CI path; Snowflake as a deployable runtime target.

### Out of scope (explicitly)

- Full-database ingestion (24M+ activities). The scoped slice targets a **low five-figure** activity count to keep build/CI fast.
- Cheminformatics modeling (QSAR/ML prediction). This is a *data-engineering* artifact; predictive modeling is a possible future extension, not part of v1.
- Bulk FTP-dump ingestion of the full SQLite/PostgreSQL release. The API-based scoped extract is the v1 path; a bulk-load path is an optional later milestone.
- Real-time / streaming ingestion. Batch, version-pinned snapshots only.

### Volume target

Scoped extract should resolve to roughly **10k–40k activities** across a few hundred to low-thousands of compounds. Exact numbers depend on the configured target set and are recorded at build time.

---

## 4. Data source — ChEMBL

- **Provider:** EMBL-EBI (European Bioinformatics Institute).
- **Release pinning:** Pin to a specific release (e.g. ChEMBL 36) and record the version in raw-layer metadata. Never silently float to "latest."
- **Access method (v1):** REST API via the official `chembl_webresource_client` Python package, with explicit pagination and retry/backoff.
- **Core entities used:**
  - `molecule` — structures (canonical SMILES, InChIKey), calculated properties, `max_phase`, `molecule_type`, `first_approval`.
  - `target` — `target_chembl_id`, `pref_name`, `target_type`, `organism`; gene/UniProt mapping via target components.
  - `assay` — `assay_chembl_id`, `assay_type`, `confidence_score`, linked target and document.
  - `activity` — `standard_type`, `standard_relation`, `standard_value`, `standard_units`, `pchembl_value`, validity flags.
- **Standardization:** ChEMBL normalizes activity endpoints and units where possible and provides `pchembl_value` (a −log10-scaled potency) for cross-assay comparability. Confidence scores annotate target–assay mapping reliability.
- **Licensing:** ChEMBL data is released under a Creative Commons Attribution-ShareAlike license (verify the exact version on the current ChEMBL site and reproduce the attribution in `README.md` and a `LICENSE-DATA` note). Provenance and version are documented in the README; no bulk data is committed to the repo beyond small fixtures.

---

## 5. Architecture

```
            ┌─────────────────────────────┐
            │   ChEMBL REST API (pinned)  │
            └──────────────┬──────────────┘
                           │  extract/  (paginated, idempotent, retry/backoff)
                           ▼
        ┌──────────────────────────────────────────┐
        │  raw schema  — landed 1:1 (Parquet)       │
        │  + ingestion metadata:                    │
        │    fetch_ts · source_endpoint ·           │
        │    chembl_version · row_hash              │   ← provenance
        └──────────────────┬───────────────────────┘
                           │  dbt staging
                           ▼
        ┌──────────────────────────────────────────┐
        │  staging schema (stg_*) — typed, cleaned, │
        │  unit-normalized, filtered, deduped       │
        └──────────────────┬───────────────────────┘
                           │  dbt marts
                           ▼
        ┌──────────────────────────────────────────┐
        │  marts schema — star schema:              │
        │  dim_compound · dim_target · dim_assay ·  │
        │  fact_activity                            │
        └──────────────────┬───────────────────────┘
                           │  dbt analytics
                           ▼
        ┌──────────────────────────────────────────┐
        │  analytics schema:                        │
        │  mart_target_sar · mart_compound_         │
        │  selectivity · mart_chemical_space        │
        └──────────────────┬───────────────────────┘
                           │
                           ▼
                 Streamlit dashboard (DuckDB-backed)
                 [optional: Tableau Public]

  Runtime engines: DuckDB (local/CI default) ⇄ Snowflake (cloud deploy)
  Orchestrated end-to-end · containerized · CI-gated
```

**Pipeline sequence (strict):** `extract → load_raw → dbt build (staging → marts → analytics) → tests → serve`. The build is idempotent and re-runnable from a pinned release.

---

## 6. Data model

### 6.0 Target-set configuration

`config/target_set.yml` is the single source of truth for the compound universe. Example shape:

```yaml
chembl_version: "36"
organism: "Homo sapiens"
targets:
  - chembl_id: CHEMBL2095173   # Tubulin (example placeholder — verify IDs)
    label: tubulin
  - chembl_id: CHEMBLxxxxxxx
    label: topoisomerase_ii
activity:
  standard_types: [IC50, GI50, Ki, Kd, EC50]
  require_pchembl: true
  min_confidence_score: 8
```

The extract layer reads this config; changing scope = changing one file (reproducible, reviewable, no code edits).

### 6.1 Raw layer (`raw`)

Landed as-is, one Parquet dataset per entity, plus ingestion metadata columns (`_fetch_ts`, `_source_endpoint`, `_chembl_version`, `_row_hash`). No business logic.

- `raw_molecules`, `raw_targets`, `raw_assays`, `raw_activities`.

### 6.2 Staging layer (`stg_*`)

| Model | Grain | Key transforms |
|---|---|---|
| `stg_molecules` | 1 row / molecule | Type properties; extract canonical SMILES + InChIKey; derive `is_approved_drug` from `max_phase`; null-handle |
| `stg_targets` | 1 row / target | Organism filter; map gene/UniProt where available; standardize `target_type` |
| `stg_assays` | 1 row / assay | Filter `confidence_score >= min_confidence`; standardize `assay_type` |
| `stg_activities` | 1 row / activity | Keep `pchembl_value not null`; `standard_relation = '='`; whitelist `standard_type`; drop rows with `data_validity_comment`; dedupe |

### 6.3 Dimensional marts (`marts`) — star schema

**`dim_compound`** (grain: 1 row / molecule)

| Column | Type | Notes |
|---|---|---|
| compound_key | int (surrogate) | PK |
| chembl_id | varchar | natural key |
| pref_name | varchar | |
| canonical_smiles | varchar | |
| inchi_key | varchar | |
| mw_freebase | double | |
| alogp | double | |
| hba / hbd | int | |
| psa | double | |
| rotatable_bonds | int | |
| num_ro5_violations | int | |
| ro3_pass | boolean | |
| max_phase | int | |
| is_approved_drug | boolean | derived |
| molecule_type | varchar | |

**`dim_target`** (grain: 1 row / target): `target_key`, `chembl_id`, `pref_name`, `target_type`, `organism`, `gene_symbol`, `uniprot_accession`.

**`dim_assay`** (grain: 1 row / assay): `assay_key`, `chembl_id`, `description`, `assay_type`, `confidence_score`.

**`fact_activity`** (grain: 1 row / measured activity):

| Column | Type | Notes |
|---|---|---|
| activity_id | bigint | natural key from ChEMBL |
| compound_key | int | FK → dim_compound |
| target_key | int | FK → dim_target |
| assay_key | int | FK → dim_assay |
| standard_type | varchar | IC50 / Ki / … |
| standard_relation | varchar | '=' after staging filter |
| standard_value | double | |
| standard_units | varchar | normalized |
| pchembl_value | double | standardized potency |
| document_chembl_id | varchar | source reference |

### 6.4 Analytical marts (`analytics`)

**`mart_target_sar`** (grain: compound × target)
Aggregates `fact_activity` per compound–target pair: `median_pchembl`, `max_pchembl`, `n_measurements`, `n_assays`. Enables potency ranking and SAR exploration per target.

**`mart_compound_selectivity`** (grain: compound)
For compounds with ≥2 targets: `n_targets`, `best_pchembl`, `best_target`, `pchembl_spread`, `selectivity_index` (best vs. second-best target potency). Surfaces selective vs. promiscuous chemotypes.

**`mart_chemical_space`** (grain: compound)
Physicochemical profile joined with potency summary and `is_approved_drug` flag, for chemical-space distribution and approved-vs-research comparison in the dashboard.

---

## 7. Tech stack & rationale

| Layer | Tool | Why / CV signal |
|---|---|---|
| Extract/Load | Python + `chembl_webresource_client` | Real API engineering: pagination, retry, idempotency |
| Validation | Pandera (ingestion) + dbt tests (warehouse) | Two-tier data quality |
| Transform | dbt | Industry-standard ELT; consistency with `mAbVault` shows depth |
| Engine (local) | DuckDB | Fast, zero-infra CI default |
| Engine (cloud) | Snowflake (free trial) | Runtime-portability target: same models, second `dbt-snowflake` profile. **Profile defined; not yet run against a live account** (see Section 11). |
| Orchestration | **Dagster** | **Implemented** as an asset graph (extract multi-asset + dbt models as assets, dbt tests as asset checks). See `docs/ORCHESTRATION.md`. |
| Serving | Streamlit (DuckDB-backed) | Functionally fills the Tableau gap; one optional Tableau Public board for the literal keyword |
| Packaging | Docker + docker-compose | `docker compose up` reproducibility |
| CI | GitHub Actions | Automated build + test gate |

**dbt adapters:** `dbt-duckdb` for local/CI, `dbt-snowflake` for the cloud target — same models, swapped profile. This is the headline "runtime-agnostic warehouse" demonstration.

---

## 8. Repository structure

```
SARVault/
├── README.md
├── SPEC.md
├── LICENSE                      # code license (e.g. MIT)
├── LICENSE-DATA.md              # ChEMBL CC BY-SA attribution
├── pyproject.toml
├── Dockerfile
├── docker-compose.yml
├── .github/workflows/ci.yml
├── config/
│   └── target_set.yml
├── extract/
│   ├── chembl_client.py         # session, retry/backoff, pagination
│   ├── extract_molecules.py
│   ├── extract_targets.py
│   ├── extract_assays.py
│   ├── extract_activities.py
│   └── load_raw.py              # write Parquet + metadata to raw
├── validation/
│   └── schemas.py               # Pandera schemas for raw entities
├── dbt/
│   ├── dbt_project.yml
│   ├── profiles/                # duckdb + snowflake profiles
│   ├── models/
│   │   ├── staging/
│   │   ├── marts/
│   │   └── analytics/
│   └── tests/
├── orchestration/
│   └── definitions.py           # Dagster (or dags/ for Airflow)
├── dashboard/
│   └── app.py                   # Streamlit
└── tests/
    ├── test_extract.py
    └── test_transforms.py
```

---

## 9. Data quality & testing

- **Ingestion (Pandera):** schema, dtype, non-null, and range checks on raw entities before landing.
- **Warehouse (dbt tests):**
  - `unique` + `not_null` on all primary/surrogate keys.
  - `relationships` (referential integrity) on every fact→dim FK.
  - `accepted_values` on `standard_type`, `target_type`.
  - Custom test: `pchembl_value` within a plausible range (e.g. 0–14).
  - Custom test: `fact_activity` grain is unique on `activity_id`.
- **Unit tests (pytest):** extract pagination logic, metadata stamping, and key transformation helpers.
- **CI gate:** `extract (cached fixture) → load_raw → dbt build → dbt test → pytest`. CI runs against DuckDB with a small committed fixture slice, not live API calls, for determinism.

---

## 10. Orchestration

A single job materializes the full lineage in the strict sequence: `extract → load_raw → dbt build → dbt test`. Assets are declared so the orchestrator's lineage view mirrors the medallion layers. Scheduling is nominal (the data is a pinned snapshot); the orchestration exists to demonstrate DAG modeling and observability, not because the source updates frequently.

---

## 11. Open decisions (to confirm before build)

1. **Orchestrator:** ~~Dagster vs. Airflow~~ → **resolved: Dagster**, implemented as an asset graph (`orchestration/definitions.py`, `docs/ORCHESTRATION.md`). Airflow remains a documented alternative.
2. **Exact target set:** Final list of ChEMBL target IDs for `target_set.yml`. Needs a short verification pass against ChEMBL to confirm IDs and resulting volume.
3. **Snowflake timing:** Trial accounts expire; sequence the Snowflake milestone so the trial window covers the demo/screenshot phase.
4. **Repo name:** ~~pending~~ → **resolved: `SARVault`** (parallels `mAbVault`; SAR = the analytical core).
5. **Bulk-load path:** Include the optional FTP-dump milestone in v1 or defer? → *Recommendation: defer; API path is sufficient for the narrative.*

---

## 12. Milestone plan (PR-based)

| Milestone | Deliverable | Branch |
|---|---|---|
| M0 | Repo scaffold, Docker, CI skeleton, `target_set.yml` | `chore/scaffold` |
| M1 | Extract layer + raw load with provenance + Pandera validation | `feat/extract-raw` |
| M2 | dbt staging models + tests | `feat/staging` |
| M3 | Dimensional marts (dims + fact) + tests | `feat/marts` |
| M4 | Analytical marts (SAR, selectivity, chemical space) | `feat/analytics` |
| M5 | Streamlit serving layer | `feat/dashboard` |
| M6 | Snowflake deploy target + dual-profile docs | `feat/snowflake` |
| M7 | Orchestration wiring | `feat/orchestration` |
| M8 *(opt)* | Tableau Public board + README polish | `docs/showcase` |

Workflow mirrors `mAbVault`: design-first, patch/PR-based, full pipeline run before commit, English code/commits.

---

## 13. Definition of done (v1)

- `docker compose up` reproduces the full warehouse from a pinned ChEMBL release.
- All dbt tests and pytest pass in CI.
- Star schema + three analytical marts materialize with documented grain.
- Dashboard renders SAR ranking, selectivity, and chemical-space views.
- README documents data provenance, ChEMBL version, license/attribution, and the DuckDB↔Snowflake runtime story.
- Snowflake build verified at least once (screenshot/notes in `docs/`). **Outstanding: never run against a live Snowflake account; the profile is defined but unverified.**
