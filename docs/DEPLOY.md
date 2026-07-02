# Deploying the SARVault dashboard (Streamlit Community Cloud)

This deploys the read-only Streamlit dashboard to a public demo URL. It does
**not** deploy the extract/dbt pipeline — the app serves a prebuilt DuckDB
warehouse. The Snowflake runtime target is a separate milestone (M6).

## How the warehouse reaches the cloud

`warehouse.duckdb` is git-ignored (no bulk ChEMBL data is committed). On a fresh
Community Cloud container the file is therefore absent, so the app fetches it
once at boot from a stable URL and caches it in the container filesystem:

- `dashboard/data.ensure_warehouse()` downloads the file only if it is missing
  **and** a URL is provided (via the `SARVAULT_WAREHOUSE_URL` secret/env). For a
  private repo it also needs `SARVAULT_WAREHOUSE_TOKEN` and fetches via the
  GitHub API; for a public repo the URL is fetched directly.
- Local runs that already have a `warehouse.duckdb` are untouched — no URL, no
  download, existing behaviour.

The recommended host for the file is a **GitHub Release asset** on this repo.

## One-time setup

### 1. Build the scoped warehouse locally

```bash
python -m extract.run
dbt build --project-dir dbt --profiles-dir dbt/profiles   # writes warehouse.duckdb
```

### 2. Publish the warehouse as a Release asset

```bash
gh release create warehouse-v1 warehouse.duckdb \
  --title "SARVault warehouse (ChEMBL 36 slice)" \
  --notes "Prebuilt DuckDB warehouse for the Streamlit demo. Rebuild via extract + dbt."
```

Copy the asset's browser download URL (Releases page → right-click the
`warehouse.duckdb` asset → copy link). It looks like:

```
https://github.com/Knusperftw/SARVault/releases/download/warehouse-v1/warehouse.duckdb
```

### 3. Deploy on Community Cloud

1. https://share.streamlit.io → **New app** → pick this repo.
2. **Branch:** `main` · **Main file path:** `dashboard/app.py`.
3. **Advanced settings → Python version:** `3.12` (matches local dev; `3.11` also
   works). Python cannot be changed after deploy without redeploying.
4. **Advanced settings → Secrets:** add

   ```toml
   SARVAULT_WAREHOUSE_URL = "https://github.com/Knusperftw/SARVault/releases/download/warehouse-v1/warehouse.duckdb"
   # Private repo only: a token with read access to this repo's contents.
   SARVAULT_WAREHOUSE_TOKEN = "github_pat_xxx"
   ```

   For a **private** repo the browser download URL 404s for anonymous requests,
   so `SARVAULT_WAREHOUSE_TOKEN` is required: the app then fetches the asset
   through the GitHub API (`GET /repos/{owner}/{repo}/releases/tags/{tag}` →
   asset id → `Accept: application/octet-stream`). Use a **fine-grained PAT**
   scoped to just this repository with **Contents: Read-only** (or a classic
   token with the `repo` scope). The token lives only in Streamlit's encrypted
   secrets, never in git. For a **public** repo, omit the token — the URL is
   fetched directly.

5. **Deploy.** First boot installs deps and downloads the warehouse (a few
   minutes); subsequent loads are fast.

## Files that make this work

| File | Purpose |
|---|---|
| `requirements.txt` | pip deps Community Cloud installs (streamlit, plotly, duckdb, rdkit, pandas, pyyaml). |
| `packages.txt` | apt deps for RDKit's 2D drawing (`libxrender1`, `libxext6`) — without these, `import rdkit ... Draw` fails with `libXrender.so.1: cannot open shared object file`. |
| `.streamlit/config.toml` | dark theme + primary colour (already present, must stay at repo root). |

## Refreshing the demo data

After a real-data refresh, rebuild and replace the asset:

```bash
python -m extract.run
dbt build --project-dir dbt --profiles-dir dbt/profiles
gh release upload warehouse-v1 warehouse.duckdb --clobber
```

Reboot the app from the Community Cloud dashboard to pull the new file (the
container caches the previous download until it restarts).

## Alternative: Git LFS (simpler, but commits the artifact)

Instead of fetch-at-boot you can track the warehouse with Git LFS — Community
Cloud pulls LFS files automatically, and no app code or secret is needed:

```bash
git lfs install
git lfs track "warehouse.duckdb"
# remove the "*.duckdb" ignore line for this file, then commit warehouse.duckdb + .gitattributes
```

Trade-off: this commits a bulk-data artifact (counter to the repo's "no bulk
ChEMBL data committed" stance) and consumes the free LFS bandwidth quota. The
fetch-at-boot path above avoids both.
