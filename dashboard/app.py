"""SARVault dashboard: sidebar navigation over the dbt-modelled warehouse.

Run from the repo root (after building the warehouse with dbt):
    streamlit run dashboard/app.py
"""

import streamlit as st

from dashboard import data
from dashboard.views import chemical_space, data_quality, overview, sar, selectivity

st.set_page_config(page_title="SARVault", layout="wide")


@st.cache_resource
def _connection():
    return data.connect()


def _scope():
    return st.session_state.get("scope")


def _overview_page():
    overview.render(_connection(), _scope())


def _sar_page():
    sar.render(_connection(), _scope())


def _selectivity_page():
    selectivity.render(_connection(), _scope())


def _chemical_space_page():
    chemical_space.render(_connection(), _scope())


def _data_quality_page():
    data_quality.render(_connection(), _scope())


# --- header (shown on every page) ---
st.markdown("# SAR:green[Vault]")
st.caption("Interactive, read-only view over the dbt-modelled ChEMBL bioactivity warehouse")

try:
    con = _connection()
except Exception as exc:  # warehouse missing or unbuilt
    st.error(
        "Could not open the warehouse. Build it first from the repo root:\n\n"
        "`dbt build --project-dir dbt --profiles-dir dbt/profiles`\n\n"
        f"Details: {exc}"
    )
    st.stop()

nav = st.navigation(
    [
        st.Page(_overview_page, title="Overview", icon="🏠", default=True),
        st.Page(_sar_page, title="SAR Ranking", icon="📊"),
        st.Page(_selectivity_page, title="Selectivity", icon="🎯"),
        st.Page(_chemical_space_page, title="Chemical Space", icon="🧪"),
        st.Page(_data_quality_page, title="Data Quality", icon="🔎"),
    ]
)

# --- global scope filter (sidebar) ---
target_names = data.list_target_names(con)
st.sidebar.divider()
st.sidebar.markdown("**🔬 Scope — targets**")
chosen = st.sidebar.multiselect(
    "Targets",
    target_names,
    default=st.session_state.get("scope", target_names),
    label_visibility="collapsed",
)
st.session_state["scope"] = chosen or target_names

nav.run()
