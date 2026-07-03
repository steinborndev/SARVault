# SPEC — ChEMBL Bioactivity SAR Warehouse

**Project:** `SARVault`
**Status:** Implemented — this document describes the **as-built** warehouse (v2 "SAR platform"). Where something is defined but not yet exercised against a live system it is called out explicitly (see §12).
**Author:** Dr. Benjamin Steinborn
**Purpose:** Portfolio project demonstrating modern data-engineering competency on a real external life-science API, complementing the internal-R&D-focused `mAbVault` project.

---

## 1. Purpose & positioning

This project builds a reproducible, layered data warehouse over the public **ChEMBL** bioactivity database, scoped to **cytotoxic / tubulin-targeting compounds** — the chemical space underlying ADC payloads and classical chemotherapeutics.

It evidences:

- Ingestion from a real external REST API (pagination, idempotency, provenance) — not flat-file reads — plus cross-reference and co-crystal-structure enrichment (UniChem, PDBe).
- A medallion-style transformation layer (raw → staging → marts) using dbt.
- A clean dimensional model with documented grain, and change-over-time modeling (an incremental fact and a Type-2 snapshot across ChEMBL releases).
- **Structural** SAR analytics on RDKit (fingerprints, similarity, scaffolds, activity cliffs, chemical-space embedding) — the domain-credible differentiator, not just potency aggregation.
- A cloud-warehouse runtime target (Snowflake) alongside a local engine (DuckDB).
- Orchestration (Dagster asset graph), data observability (source freshness, anomaly guards, published dbt docs), containerization, and CI.

### Relationship to `mAbVault`

| Axis | `mAbVault` | This project |
|---|---|---|
| Data origin | Simulated internal R&D | Real external public API |
| Domain | mAb stability / CMC | Cytotoxic small-molecule SAR |
| Warehouse | DuckDB (local) | DuckDB **and** Snowflake-ready (cloud) |
| Emphasis | Data *content* understanding | Modern DE *stack* breadth + structural SAR |

Together they signal: *"I understand pharma R&D data and I can build the modern stack against real sources."*

### Target roles

R&D Data Engineering (Roche) and CMC / analytical data roles (e.g. Tubulis). The payload-chemistry scope deliberately ties to the author's B.Sc. thesis (chondramides / microtubule targeting), Ph.D. (chemotherapeutic delivery), and ADC-payload relevance.

---

## 2. Analytical thesis

A pure "copy ChEMBL into a warehouse" pipeline has no narrative. This warehouse answers concrete questions about the payload chemical space, in two tiers.

**Potency & selectivity (per target):**

1. **SAR ranking** — For a given target (tubulin / β-tubulin, topoisomerase I / II-α), which compounds are most potent, and how is potency distributed?
2. **Selectivity profiling** — For compounds tested against multiple targets, what is the potency spread and selectivity index?
3. **Chemical-space characterization** — How are physicochemical properties (MW, logP, HBA/HBD, Ro5) distributed, and where do approved drugs sit?

**Structural SAR (per chemotype):**

4. **Structural analogs** — Given a compound, its nearest neighbours by ECFP4 Tanimoto, annotated with the potency delta (a lead signal).
5. **Chemical series** — Compounds grouped by Bemis-Murcko scaffold, with per-series size, potency spread and target reach.
6. **Activity cliffs** — Pairs of highly similar compounds whose potency differs sharply on the same target (ranked by SALI) — the sharpest signal in SAR.
7. **Structural embedding** — A 2-D UMAP of ECFP4 fingerprints showing where potency concentrates in chemical space.

These questions are materialized as the analytical marts in §6.4.

---

## 3. Scope

### In scope

- A **config-driven target/mechanism set** (§6.0) defining the compound universe: tubulin, β-tubulin, topoisomerase I and II-α.
- ChEMBL core entities: molecules, targets, assays, activities. Enrichment: UniChem cross-references, PDBe co-crystal structures.
- Human-target focus by default (`organism = Homo sapiens`), configurable.
- Activity measurements with a standardized potency value (`pchembl_value`).
- **Deterministic structural cheminformatics** on RDKit: ECFP4 fingerprints, Bemis-Murcko scaffolds, Tanimoto similarity, substructure search, activity cliffs, UMAP embedding.
- Local (DuckDB) build as the primary CI path; Snowflake as a defined, deployable runtime target.

### Out of scope (explicitly)

- **QSAR / ML potency prediction.** This is a *data-engineering* artifact; the structural analytics above are deterministic, not learned. Predictive modeling is a possible v3 extension, not part of this scope.
- Full-database ingestion (24M+ activities). The scoped slice keeps build/CI fast.
- Bulk FTP-dump ingestion of the full ChEMBL SQLite/PostgreSQL release. (UniChem *is* consumed via its bulk FTP dumps; the ChEMBL core is the scoped REST extract.)
- Real-time / streaming ingestion. Batch, version-pinned snapshots only.

### Volume target

The scoped extract resolves to a **low five-figure** activity count across a few hundred to low-thousands of compounds. Exact numbers are recorded at build time (`scripts/profile_sar.py` → `docs/DATA_PROFILE.md`).

---

## 4. Data sources

### ChEMBL (primary)

- **Provider:** EMBL-EBI. **Release pinning:** pinned via `config/target_set.yml` (`chembl_version`, currently **ChEMBL 36**) and stamped into raw-layer metadata; never silently floats to "latest."
- **Access:** the ChEMBL **REST API** via a custom `requests` session (`extract/chembl_client.py`) with explicit pagination (offset, chunk size ≤ 50 on ID filters) and retry/backoff — not the `chembl_webresource_client` wrapper.
- **Core entities:** `molecule` (SMILES, InChIKey, calculated properties, `max_phase`, `molecule_type`), `target` (`target_chembl_id`, `pref_name`, `target_type`, `organism`), `assay` (`assay_chembl_id`, `assay_type`, `confidence_score`), `activity` (`standard_type/relation/value/units`, `pchembl_value`, validity flags).
- **Standardization:** ChEMBL provides `pchembl_value` (−log10-scaled potency) for cross-assay comparability; confidence scores annotate target–assay mapping reliability.
- **Licensing:** CC BY-SA (see `LICENSE-DATA.md`). Provenance and version are documented in the README; no bulk data is committed beyond small fixtures.

### Enrichment

- **UniChem** — compound cross-references (PubChem, DrugBank, PDBe, BindingDB, SureChEMBL) from UniChem's **bulk FTP dumps** (the live API returned 5xx; the dumps are the reliable path).
- **PDBe** — co-crystal structures per ligand via per-id REST GET requests (`in_pdb` resolution + entry summary: title, method, year, resolution) feeding the embedded 3-D Mol\* viewer.

### Derived (in-pipeline, RDKit)

- **Cheminfo stage** (`extract/cheminfo.py`) — for every parseable SMILES, an ECFP4 (Morgan radius 2, 2048-bit) fingerprint and a Bemis-Murcko scaffold, landed as raw Parquet with provenance (`_source = rdkit`, pinned `_rdkit_version`). Deterministic given a pinned RDKit.
- **Embedding stage** (`extract/embedding.py`) — a fixed-seed 2-D UMAP projection of the ECFP4 fingerprints, precomputed once so the dashboard only plots.

---

## 5. Architecture

```
   ChEMBL REST API (pinned)      UniChem bulk FTP      PDBe REST
            │                          │                  │
            │  extract/ (paginated, idempotent, retry/backoff)
            ▼                          ▼                  ▼
   ┌─────────────────────────────────────────────────────────────┐
   │  raw schema — landed as Parquet + provenance                │
   │  (_fetch_ts · _source · _chembl_version · _row_hash)        │
   │  core: molecules · targets · assays · activities            │
   │  enrich: xref_unichem · pdbe_structures · pdbe_summary       │
   │  derived (RDKit): compound_cheminfo · compound_embedding    │
   └──────────────────────────────┬──────────────────────────────┘
                                  │  dbt staging (materialized: table)
                                  ▼
   ┌─────────────────────────────────────────────────────────────┐
   │  staging (stg_*) — typed, cleaned, filtered, deduped        │
   └──────────────────────────────┬──────────────────────────────┘
                                  │  dbt marts
                                  ▼
   ┌─────────────────────────────────────────────────────────────┐
   │  marts — star schema:                                       │
   │  dim_compound · dim_target · dim_assay · dim_scaffold       │
   │  fact_activity (incremental)   + snapshot: compound_status  │
   │                                  (SCD2 across ChEMBL rels.)  │
   └──────────────────────────────┬──────────────────────────────┘
                                  │  dbt analytics
                                  ▼
   ┌─────────────────────────────────────────────────────────────┐
   │  analytics — 9 marts: target_sar · compound_selectivity ·   │
   │  chemical_space · compound_fingerprint · chemical_series ·  │
   │  activity_cliff · compound_pdb · compound_xref ·            │
   │  compound_catalog                                           │
   └──────────────────────────────┬──────────────────────────────┘
                                  ▼
                   Streamlit dashboard (DuckDB-backed)

  Orchestration:  Dagster asset graph (extract → raw → dbt build → tests),
                  dbt tests + explicit checks surfaced as asset checks.
  Observability:  dbt source freshness · anomaly guards · published dbt docs.
  Runtime:        DuckDB (local/CI default)  ⇄  Snowflake-ready (second profile).
  Delivery:       containerized · CI-gated · warehouse shipped via GitHub Release.
```

**Pipeline sequence (strict):** `extract → load_raw → cheminfo/embedding → dbt build (staging → marts → analytics + snapshot) → tests → serve`. The build is idempotent and re-runnable from a pinned release; `fact_activity` refreshes incrementally.

---

## 6. Data model

### 6.0 Target-set configuration

`config/target_set.yml` is the single source of truth for the compound universe. Actual shape:

```yaml
chembl_version: "36"
organism: "Homo sapiens"
targets:
  - chembl_id: "CHEMBL2095182"   # Tubulin (PROTEIN COMPLEX GROUP)
    label: "tubulin"
  - chembl_id: "CHEMBL3832942"   # Tubulin beta (PROTEIN FAMILY)
    label: "tubulin_beta"
  - chembl_id: "CHEMBL1781"      # DNA topoisomerase 1 (TOP1)
    label: "topoisomerase_i"
  - chembl_id: "CHEMBL1806"      # DNA topoisomerase 2-alpha (TOP2A)
    label: "topoisomerase_ii_alpha"
activity:
  standard_types: [IC50, GI50, Ki, Kd, EC50]
  require_pchembl: true
  min_confidence_score: 5   # permissive: the tubulin COMPLEX/FAMILY targets carry
                            # lower confidence than the single-protein topos
```

Changing scope = editing this file only (reproducible, reviewable, no code edits).

### 6.1 Raw layer (`raw`)

Landed as-is, one Parquet dataset per entity, plus provenance columns (`_fetch_ts`, `_source`, `_chembl_version` / `_rdkit_version`, `_row_hash`). No business logic.

- Core: `raw_molecules`, `raw_targets`, `raw_assays`, `raw_activities`.
- Enrichment: `raw_xref_unichem`, `raw_pdbe_structures`, `raw_pdbe_summary`.
- Derived (RDKit): `raw_compound_cheminfo`, `raw_compound_embedding`.

### 6.2 Staging layer (`stg_*`, materialized as tables)

| Model | Grain | Key transforms |
|---|---|---|
| `stg_molecules` | 1 / molecule | Type properties; canonical SMILES + InChIKey; derive `is_approved_drug` from `max_phase` |
| `stg_targets` | 1 / target | Organism filter; standardize `target_type` |
| `stg_assays` | 1 / assay | Filter `confidence_score >= min_confidence`; standardize `assay_type` |
| `stg_activities` | 1 / activity | `pchembl_value not null`; `standard_relation = '='`; whitelist `standard_type`; drop `data_validity_comment`; dedupe |
| `stg_compound_cheminfo` | 1 / compound | ECFP4 fingerprint (hex + on-bit list) and Murcko scaffold |
| `stg_compound_embedding` | 1 / compound | 2-D UMAP coordinates of the ECFP4 fingerprint |
| `stg_compound_xref` | 1 / (compound, ref) | Explode ChEMBL `cross_references` JSON |
| `stg_compound_xref_unichem` | 1 / (compound, ref) | UniChem bulk cross-references, scoped to our compounds |
| `stg_pdbe_structure` | 1 / (ligand, PDB entry) | PDBe `in_pdb` resolution |
| `stg_pdbe_summary` | 1 / PDB id | PDBe entry metadata (title, method, year, resolution) |

Staging is materialized as **tables** (not views) so the shipped `warehouse.duckdb` is self-contained — a deployment ships only the `.duckdb`, without the raw Parquet a view would re-read.

### 6.3 Dimensional marts (`marts`) — star schema

**`dim_compound`** (1 / molecule): `compound_key` (surrogate PK), `molecule_chembl_id` (natural key), `pref_name`, `canonical_smiles`, `inchi_key`, `mw_freebase`, `alogp`, `hba`, `hbd`, `psa`, `rotatable_bonds`, `num_ro5_violations`, `ro3_pass`, `aromatic_rings`, `qed_weighted`, `max_phase`, `is_approved_drug`, `molecule_type`.

**`dim_target`** (1 / target): `target_key`, `target_chembl_id`, `pref_name`, `target_type`, `organism`, …

**`dim_assay`** (1 / assay): `assay_key`, `assay_chembl_id`, `description`, `assay_type`, `confidence_score`.

**`dim_scaffold`** (1 / distinct Bemis-Murcko scaffold): `scaffold_key`, `murcko_scaffold_smiles`. The structural grouping dimension behind chemical series.

**`fact_activity`** (1 / measured activity) — **incremental**, `unique_key = activity_id`, `delete+insert`. Columns: `activity_id`, `compound_key`/`target_key`/`assay_key` (FKs), `standard_type`, `standard_relation`, `standard_value`, `standard_units`, `pchembl_value`, `document_chembl_id`. Activities are append-only with a stable id, so a refresh loads only the delta; a no-delta re-run is a no-op.

**`compound_status`** (dbt **snapshot**, SCD2) — Type-2 history of a compound's development status (`max_phase`, `is_approved_drug`, `pref_name`), keyed on the ChEMBL natural key `molecule_chembl_id` (not the surrogate `compound_key`, which is a non-release-stable `row_number()`). A new validity window opens when a payload's approval status changes across ChEMBL releases.

### 6.4 Analytical marts (`analytics`)

| Mart | Grain | Purpose |
|---|---|---|
| `mart_target_sar` | compound × target | `median/max_pchembl`, `n_measurements`, `n_assays` — potency ranking per target |
| `mart_compound_selectivity` | compound (≥2 targets) | `n_targets`, `best_pchembl`, `pchembl_spread`, `selectivity_index` |
| `mart_chemical_space` | compound | Physicochemical profile + potency + UMAP coordinates + approval flag |
| `mart_compound_fingerprint` | compound | ECFP4 fingerprint + scaffold key — the substrate for similarity/series/cliffs |
| `mart_chemical_series` | scaffold | Per-series size, potency spread, target reach |
| `mart_activity_cliff` | (target, compound_a, compound_b) | Similar pairs with sharp Δpotency, ranked by SALI |
| `mart_compound_pdb` | compound | Co-crystal PDB availability + metadata for the 3-D viewer |
| `mart_compound_xref` | compound | Consolidated external cross-references |
| `mart_compound_catalog` | compound | Library-page catalog joining structure, potency, Ro5, PDB flags |

Cliff thresholds are dbt vars (`cliff_min_tanimoto`, `cliff_min_delta_pchembl`), chosen from the real-data profile.

---

## 7. Tech stack & rationale

| Layer | Tool | Why / CV signal |
|---|---|---|
| Extract/Load | Python + `requests` (custom session) | Real API engineering: pagination, retry, idempotency, provenance |
| Enrichment | UniChem (bulk FTP), PDBe (REST) | Cross-reference + co-crystal-structure integration |
| Cheminformatics | **RDKit** (+ UMAP) | ECFP4 fingerprints, Murcko scaffolds, Tanimoto, substructure, activity cliffs, 2-D embedding |
| Validation | Pandera (ingestion) + dbt tests (warehouse) | Two-tier data quality |
| Transform | dbt (medallion) | Industry-standard ELT; incremental + SCD2 snapshot for change-over-time |
| Engine (local) | DuckDB | Fast, zero-infra CI default |
| Engine (cloud) | Snowflake | Runtime-portability target: same models, second `dbt-snowflake` profile. **Profile defined; not yet run against a live account** (see §12). |
| Orchestration | **Dagster** (`dagster-dbt`) | Asset graph mirroring the medallion; dbt tests + explicit checks as asset checks (`docs/ORCHESTRATION.md`) |
| Observability | dbt source freshness · anomaly guards · published dbt docs | "Is the data right today?" — freshness, row-count/null-rate/distribution guards, lineage site |
| Serving | Streamlit (DuckDB-backed) | Multipage dashboard, deployed on Streamlit Community Cloud; embedded Mol\* 3-D viewer |
| Packaging | Docker + docker-compose | `docker compose up` reproducibility |
| CI | GitHub Actions | ruff → dbt build → source freshness → dagster job → pytest; docs site to GitHub Pages |

**dbt adapters:** `dbt-duckdb` for local/CI, `dbt-snowflake` for the cloud target — same models, swapped profile. This is the headline "runtime-agnostic warehouse" demonstration (the live Snowflake proof is the one outstanding item, §12).

---

## 8. Repository structure

```
SARVault/
├── README.md · SPEC.md · LICENSE · LICENSE-DATA.md · pyproject.toml
├── Dockerfile · docker-compose.yml · requirements.txt · packages.txt
├── .github/workflows/ci.yml
├── config/target_set.yml
├── extract/
│   ├── chembl_client.py        # session, retry/backoff, pagination
│   ├── config.py · run.py
│   ├── extract_{molecules,targets,assays,activities}.py
│   ├── load_raw.py             # Parquet + provenance → raw
│   ├── unichem.py · pdbe.py    # enrichment
│   └── cheminfo.py · embedding.py   # RDKit fingerprints/scaffolds; UMAP
├── validation/schemas.py       # Pandera schemas
├── dbt/
│   ├── dbt_project.yml · profiles/   # duckdb + snowflake
│   ├── models/{staging,marts,analytics}/
│   ├── snapshots/compound_status.sql # SCD2
│   └── tests/                  # singular tests incl. anomaly guards
├── orchestration/definitions.py      # Dagster asset graph
├── dashboard/
│   ├── app.py · data.py · logic.py · chem.py · viewer.py …
│   └── views/                  # library, SAR, selectivity, series, cliffs, space, DQ
├── scripts/profile_sar.py      # real-data profiler → docs/DATA_PROFILE.md
├── docs/{ORCHESTRATION,DEPLOY,DATA_PROFILE}.md
└── tests/                      # pytest suite + fixtures/ (raw, raw_v2 release delta)
```

---

## 9. Data quality & testing

- **Ingestion (Pandera):** schema, dtype, non-null and range checks on raw entities before landing.
- **Warehouse (dbt tests):** `unique` + `not_null` on all keys; `relationships` on every fact→dim FK; `accepted_values` on `standard_type`/`target_type`; singular tests for grain (SAR, cliff, PDB, xref, snapshot), `pchembl` range, selectivity non-negativity, and series/scaffold reconciliation.
- **Change-over-time tests:** incremental no-op / delta-insert behavior and SCD2 validity-window integrity, exercised across a v1 fixture and a `raw_v2` "next release" delta (`tests/test_snapshot.py`).
- **Observability:** `dbt source freshness` over the raw layer (from `_fetch_ts`); three anomaly guards — per-layer **row-count floor**, `canonical_smiles` **null-rate ceiling**, `pchembl` **mean band** — proven to fire on a seeded bad fixture (`tests/test_observability.py`).
- **Unit tests (pytest):** extract pagination / provenance, RDKit fingerprint determinism, similarity / substructure, cliff detection, embedding reproducibility, dashboard logic, orchestration graph.
- **CI gate:** `ruff → dbt build (fixtures) → dbt source freshness → dbt parse → dagster job execute → pytest`, all against DuckDB with a committed fixture slice (no live API calls, for determinism). The Dagster job materializes the full asset graph and runs every dbt test plus explicit warehouse checks as asset checks. On `main`, CI also publishes the dbt docs site to GitHub Pages.

---

## 10. Orchestration

The full lineage is a **Dagster asset graph** (`orchestration/definitions.py`): the ChEMBL extract is a multi-asset publishing the raw roots, each dbt model is an asset via `dagster-dbt`, and dbt tests plus two explicit warehouse checks (non-empty `fact_activity`, `pchembl ∈ [0,14]`) surface as asset checks. `dagster dev` renders the lineage graph; `dagster job execute -j sarvault_transform` materializes the warehouse end-to-end from fixtures in CI. Scheduling is nominal (the source is a pinned snapshot); the orchestration exists to demonstrate DAG modeling and observability. See `docs/ORCHESTRATION.md`.

---

## 11. Design decisions (resolved)

1. **Orchestrator:** Dagster (over Airflow) — implemented as an asset graph. Airflow remains a documented alternative.
2. **Repo name:** `SARVault` (parallels `mAbVault`; SAR = the analytical core).
3. **Extract path:** scoped REST extract via `requests` over the `chembl_webresource_client` wrapper, for explicit control of pagination, provenance and retry.
4. **Bulk-load path:** ChEMBL core stays on the scoped REST extract; only UniChem is consumed via bulk FTP (its live API was unreliable).
5. **Cheminfo compute location:** in-pipeline (RDKit stage landed as raw Parquet), not on-the-fly in Streamlit, so the dashboard stays fast and the fingerprints are a governed asset.
6. **Snapshot key:** SCD2 keyed on `molecule_chembl_id` (natural key), because the surrogate `compound_key` is a `row_number()` and not release-stable.

---

## 12. Delivery status

**v1 — foundational warehouse (complete):** config-driven ChEMBL extract with provenance; dbt medallion (staging → star-schema marts → analytical marts) on DuckDB; Pandera + dbt tests; Streamlit dashboard; Docker + CI.

**Epic 1 — SAR intelligence layer (complete):** RDKit cheminfo stage (`dim_scaffold`, `mart_compound_fingerprint`); structural analogs (Tanimoto); substructure / SMARTS search; chemical series; activity cliffs (SALI); UMAP chemical-space embedding.

**Epic 2 — production DE backbone:**
- Dagster wired for real (asset graph + asset checks) — **done**.
- Incremental `fact_activity` + SCD2 `compound_status` snapshot across ChEMBL releases — **done**.
- Observability: source freshness, anomaly guards, published dbt docs site — **done**.
- **Snowflake proof — outstanding:** the dual-runtime story is real only once the same models build green against a live Snowflake account (a secrets-gated CI run + `docs/SNOWFLAKE.md` parity notes). The profile is defined; the run is not yet done.

**Epic 3 — data-as-a-product serving layer (future):** a FastAPI read API over the marts, a metrics/semantic layer, and a versioned, checksummed data release package.

**Explicit non-goal:** QSAR / ML potency prediction — kept out to preserve the data-engineering framing (a possible v3).

---

## 13. Definition of done

**v1 (met):**

- `docker compose up` reproduces the full warehouse from a pinned ChEMBL release.
- All dbt tests and pytest pass in CI.
- Star schema + analytical marts materialize with documented grain.
- Dashboard renders SAR ranking, selectivity, and chemical-space views.
- README documents data provenance, ChEMBL version, license/attribution, and the DuckDB↔Snowflake runtime story.

**Outstanding for the "runtime-agnostic" claim:** a verified live Snowflake build (dated run + parity notes in `docs/`). The profile is defined but has not been run against a live account.
