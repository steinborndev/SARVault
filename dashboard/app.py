"""SARVault dashboard: sidebar navigation over the dbt-modelled warehouse.

Run from the repo root (after building the warehouse with dbt):
    streamlit run dashboard/app.py
"""

from pathlib import Path

import base64
import os
import sys

# When launched as `streamlit run dashboard/app.py` (e.g. on Streamlit Community
# Cloud, where this package isn't pip-installed), the entrypoint's own directory
# is on sys.path but the repo root isn't — put it there before the first-party
# imports below so `dashboard` and `extract` resolve.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import streamlit as st  # noqa: E402

from dashboard import data  # noqa: E402
from dashboard.views import (  # noqa: E402
    chemical_space,
    compound_library,
    data_quality,
    overview,
    sar,
    selectivity,
)

_ASSETS = Path(__file__).resolve().parents[1] / "assets"
_LOGO = str(_ASSETS / "logo.svg")
_ICON = str(_ASSETS / "icon.svg")
_APPROVAL = ["all", "approved", "research"]

st.set_page_config(page_title="SARVault", page_icon=_ICON, layout="wide")


def _logo_html(height: int = 54) -> str:
    b64 = base64.b64encode(Path(_LOGO).read_bytes()).decode()
    return f'<img src="data:image/svg+xml;base64,{b64}" height="{height}" style="display:block">'


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


def _chemical_space_page():
    chemical_space.render(_connection(), _scope())


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

# --- header: logo (top left) + scope (top right) ---
header_left, header_right = st.columns([12, 1])
with header_left:
    st.markdown(_logo_html(), unsafe_allow_html=True)
    st.caption(
        "Interactive, read-only view over the dbt-modelled ChEMBL bioactivity warehouse"
    )
target_names = data.list_target_names(con)
prev = st.session_state.get("scope", {})
with header_right.popover("Scope", width="stretch"):
    sel_targets = st.multiselect("Targets", target_names, default=prev.get("targets", target_names))
    approval = st.radio(
        "Approval", _APPROVAL, index=_APPROVAL.index(prev.get("approval", "all")), horizontal=True
    )
    min_pchembl = st.slider("Min best pChEMBL", 0.0, 12.0, float(prev.get("min_pchembl", 0.0)), 0.5)
    structure_only = st.toggle(
        "Only with 3D crystal structure", value=bool(prev.get("structure_only", False))
    )
st.session_state["scope"] = {
    "targets": sel_targets or target_names,
    "approval": approval,
    "min_pchembl": min_pchembl,
    "structure_only": structure_only,
}

nav = st.navigation(
    [
        st.Page(_overview_page, title="Overview", icon="🏠", default=True),
        st.Page(_library_page, title="Compound Library", icon="📚"),
        st.Page(_sar_page, title="SAR Ranking", icon="📊"),
        st.Page(_selectivity_page, title="Selectivity", icon="🎯"),
        st.Page(_chemical_space_page, title="Chemical Space", icon="🧪"),
        st.Page(_data_quality_page, title="Data Quality", icon="🔎"),
    ]
)
nav.run()
