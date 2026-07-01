"""Compound library: filterable list + per-compound detail drill-down."""

import base64
import math

import pandas as pd
import streamlit as st

from dashboard import charts, chem, data, logic

_CHEMBL_URL = "https://www.ebi.ac.uk/chembl/compound_report_card/{}/"

_LIST_COLS = [
    "molecule_chembl_id",
    "pref_name",
    "is_approved_drug",
    "best_pchembl",
    "best_target",
    "n_targets",
    "mw_freebase",
    "alogp",
    "num_ro5_violations",
]

_PROP_ROWS = [
    ("MW (freebase)", "mw_freebase"),
    ("logP (AlogP)", "alogp"),
    ("HBA", "hba"),
    ("HBD", "hbd"),
    ("PSA", "psa"),
    ("Rotatable bonds", "rotatable_bonds"),
    ("Aromatic rings", "aromatic_rings"),
    ("Ro5 violations", "num_ro5_violations"),
    ("Ro3 pass", "ro3_pass"),
    ("QED", "qed_weighted"),
    ("Max phase", "max_phase"),
]


def _fmt(value) -> str:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return "—"
    if isinstance(value, float):
        return f"{value:.2f}"
    return str(value)


def _structure_html(svg: str) -> str:
    b64 = base64.b64encode(svg.encode()).decode()
    return f'<img src="data:image/svg+xml;base64,{b64}" width="320">'


def render(con, scope):
    st.header("📚 Compound library")
    target_sar = data.load_target_sar(con)
    catalog = data.load_compound_catalog(con)
    keys = logic.resolve_scope_keys(target_sar, catalog, scope)
    cat = catalog[catalog["compound_key"].isin(keys)]
    if cat.empty:
        st.info("No compounds in the current scope.")
        return

    with st.sidebar:
        st.markdown("### Filters")
        query = st.text_input("Search (ChEMBL ID / name)").strip().lower()
        p_lo, p_hi = st.slider("best pChEMBL", 0.0, 14.0, (0.0, 14.0), step=0.5)
        mw_series = cat["mw_freebase"].dropna()
        mw_top = mw_series.max() if not mw_series.empty else 1000.0
        mw_ceiling = int(math.ceil(mw_top / 50) * 50)
        mw_lo, mw_hi = st.slider("MW range", 0, mw_ceiling, (0, mw_ceiling), step=50)
        logp_lo, logp_hi = st.slider("logP range", -5.0, 12.0, (-5.0, 12.0), step=0.5)
        max_ro5 = st.slider("Max Ro5 violations", 0, 4, 4)

    # between() is inclusive on both ends, so a compound on a boundary (e.g. MW 500)
    # matches both the 450-500 and 500-550 windows.
    view = cat[
        cat["best_pchembl"].between(p_lo, p_hi)
        & cat["mw_freebase"].fillna(0).between(mw_lo, mw_hi)
        & cat["alogp"].fillna(0).between(logp_lo, logp_hi)
        & (cat["num_ro5_violations"].fillna(0) <= max_ro5)
    ]
    if query:
        view = view[
            view["molecule_chembl_id"].str.lower().str.contains(query)
            | view["pref_name"].fillna("").str.lower().str.contains(query)
        ]

    disp = view.sort_values("best_pchembl", ascending=False).reset_index(drop=True)
    st.caption(f"{len(disp)} compounds — click a row to inspect it")
    event = st.dataframe(
        disp[_LIST_COLS],
        hide_index=True,
        width="stretch",
        on_select="rerun",
        selection_mode="single-row",
        key="lib_rows",
    )
    if disp.empty:
        return

    # a clicked row drives the selection in "Inspect a compound"
    selected_rows = event.selection.rows
    if selected_rows and selected_rows[0] < len(disp):
        st.session_state["inspect_compound"] = disp.iloc[selected_rows[0]]["molecule_chembl_id"]

    st.subheader("Inspect a compound")
    options = disp["molecule_chembl_id"].tolist()
    if st.session_state.get("inspect_compound") not in options:
        st.session_state["inspect_compound"] = options[0]
    chosen = st.selectbox("Compound", options, key="inspect_compound", label_visibility="collapsed")
    row = disp[disp["molecule_chembl_id"] == chosen].iloc[0]

    st.markdown(f"### {chosen} — {row['pref_name'] or chosen}")
    st.link_button("View on ChEMBL", _CHEMBL_URL.format(chosen))

    left, right = st.columns(2)
    with left:
        svg = chem.smiles_to_svg(row.get("canonical_smiles"))
        if svg:
            st.markdown(_structure_html(svg), unsafe_allow_html=True)
        else:
            st.info("No structure available.")
        props = pd.DataFrame(
            [{"property": label, "value": _fmt(row.get(col))} for label, col in _PROP_ROWS]
        )
        st.dataframe(props, hide_index=True, width="stretch")
    with right:
        st.markdown("**Per-target potency**")
        profile = data.compound_target_profile(con, int(row["compound_key"]))
        st.plotly_chart(charts.compound_potency_bar(profile), width="stretch")
        st.dataframe(profile, hide_index=True, width="stretch")
        if row["n_targets"] >= 2 and pd.notna(row["selectivity_index"]):
            st.metric("Selectivity index (log10 fold)", round(float(row["selectivity_index"]), 2))
