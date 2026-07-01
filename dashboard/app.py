"""SARVault dashboard: sidebar navigation over the dbt-modelled warehouse.

Run from the repo root (after building the warehouse with dbt):
    streamlit run dashboard/app.py
"""

from pathlib import Path

import base64

import streamlit as st

from dashboard import data
from dashboard.views import (
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
st.session_state["scope"] = {
    "targets": sel_targets or target_names,
    "approval": approval,
    "min_pchembl": min_pchembl,
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
