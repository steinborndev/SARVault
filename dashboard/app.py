"""SARVault dashboard: sidebar navigation over the dbt-modelled warehouse.

Run from the repo root (after building the warehouse with dbt):
    streamlit run dashboard/app.py
"""

from pathlib import Path

import os
import sys

# When launched as `streamlit run dashboard/app.py` (e.g. on Streamlit Community
# Cloud, where this package isn't pip-installed), the entrypoint's own directory
# is on sys.path but the repo root isn't — put it there before the first-party
# imports below so `dashboard` and `extract` resolve.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import streamlit as st  # noqa: E402

from dashboard import branding, data, scope  # noqa: E402
from dashboard.views import (  # noqa: E402
    activity_cliffs,
    chemical_series,
    chemical_space,
    compound_library,
    data_quality,
    overview,
    sar,
    selectivity,
)

_ASSETS = Path(__file__).resolve().parents[1] / "assets"
_ICON = str(_ASSETS / "icon.svg")

st.set_page_config(page_title="SARVault", page_icon=_ICON, layout="wide")


@st.cache_resource
def _connection():
    return data.connect()


def _warehouse_url() -> str | None:
    """Warehouse download URL from Streamlit secrets or env (cloud deploy only)."""
    try:
        url = st.secrets.get("SARVAULT_WAREHOUSE_URL")
    except Exception:  # no secrets file locally
        url = None
    return url or os.environ.get("SARVAULT_WAREHOUSE_URL")


def _warehouse_token() -> str | None:
    """GitHub token for fetching the warehouse asset from a private repo."""
    try:
        tok = st.secrets.get("SARVAULT_WAREHOUSE_TOKEN")
    except Exception:
        tok = None
    return tok or os.environ.get("SARVAULT_WAREHOUSE_TOKEN")


def _scope():
    return st.session_state.get("scope", {})


def _overview_page():
    overview.render(_connection(), _scope())


def _library_page():
    compound_library.render(_connection(), _scope())


def _sar_page():
    sar.render(_connection(), _scope())


def _selectivity_page():
    selectivity.render(_connection(), _scope())


def _cliffs_page():
    activity_cliffs.render(_connection(), _scope())


def _chemical_space_page():
    chemical_space.render(_connection(), _scope())


def _chemical_series_page():
    chemical_series.render(_connection(), _scope())


def _data_quality_page():
    data_quality.render(_connection(), _scope())


try:
    data.ensure_warehouse(url=_warehouse_url(), token=_warehouse_token())
    con = _connection()
except Exception as exc:  # warehouse missing or unbuilt
    st.error(
        "Could not open the warehouse. Build it first from the repo root:\n\n"
        "`dbt build --project-dir dbt --profiles-dir dbt/profiles`\n\n"
        f"Details: {exc}"
    )
    st.stop()

# --- header: logo flush-left and the scope popover share the same top line ---
target_names = data.list_target_names(con)
logo_col, scope_col = st.columns([11, 1], vertical_alignment="top")
with logo_col:
    st.markdown(branding.logo_html(), unsafe_allow_html=True)
with scope_col:
    scope.render(target_names)

st.caption(
    "Structure–activity intelligence over a reproducible warehouse of public ChEMBL "
    "bioactivity data, scoped to cytotoxic / tubulin-targeting compounds: the chemistry "
    "behind ADC payloads and classical chemotherapeutics."
)

nav = st.navigation(
    [
        st.Page(_overview_page, title="Overview", default=True),
        st.Page(_library_page, title="Compound Library"),
        st.Page(_sar_page, title="SAR Ranking"),
        st.Page(_cliffs_page, title="Activity Cliffs"),
        st.Page(_chemical_series_page, title="Chemical Series"),
        st.Page(_selectivity_page, title="Selectivity"),
        st.Page(_chemical_space_page, title="Chemical Space"),
        st.Page(_data_quality_page, title="Data Quality"),
    ]
)
nav.run()
