"""Compound library: filterable list + per-compound detail drill-down."""

import math

import streamlit as st

from dashboard import chem, compound_detail, data, logic

_LIST_COLS = [
    "molecule_chembl_id",
    "pref_name",
    "is_approved_drug",
    "has_pdb",
    "best_pchembl",
    "best_target",
    "n_targets",
    "mw_freebase",
    "alogp",
    "num_ro5_violations",
]


def _list_column_config():
    """Display config for the compound list: formatted numbers, flags, potency bar."""
    return {
        "molecule_chembl_id": st.column_config.TextColumn("ChEMBL ID"),
        "pref_name": st.column_config.TextColumn("Name"),
        "is_approved_drug": st.column_config.CheckboxColumn("Approved"),
        "has_pdb": st.column_config.CheckboxColumn("3D", help="Has a co-crystal PDB entry"),
        "best_pchembl": st.column_config.ProgressColumn(
            "Best pChEMBL", format="%.2f", min_value=0.0, max_value=14.0
        ),
        "best_target": st.column_config.TextColumn("Best target"),
        "n_targets": st.column_config.NumberColumn("Targets", format="%d"),
        "mw_freebase": st.column_config.NumberColumn("MW", format="%.0f"),
        "alogp": st.column_config.NumberColumn("logP", format="%.1f"),
        "num_ro5_violations": st.column_config.NumberColumn("Ro5 viol.", format="%d"),
    }


def _int_max(series, default: int) -> int:
    clean = series.dropna()
    return int(clean.max()) if not clean.empty else default


@st.cache_data(show_spinner=False)
def _substructure_hits(smarts: str, smiles: tuple[str, ...]) -> dict:
    """Map each unique SMILES -> whether it contains the SMARTS (cached; RDKit is the cost)."""
    return {s: bool(chem.has_substructure(s, smarts)) for s in set(smiles) if s}


def render(con, scope):
    st.header("Compound library")
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
        hbd_max = _int_max(cat["hbd"], 10)
        hbd_lo, hbd_hi = st.slider("HBD range", 0, max(hbd_max, 1), (0, max(hbd_max, 1)))
        hba_max = _int_max(cat["hba"], 15)
        hba_lo, hba_hi = st.slider("HBA range", 0, max(hba_max, 1), (0, max(hba_max, 1)))
        max_ro5 = st.slider("Max Ro5 violations", 0, 4, 4)
        smarts_raw = st.text_input(
            "Substructure (SMARTS)",
            help="e.g. c1ccccc1 (benzene), C(=O)N (amide), [#7] (any nitrogen)",
        ).strip()

    # Validate the SMARTS once; an invalid pattern warns and is ignored (no crash).
    smarts = None
    if smarts_raw:
        if chem.is_valid_smarts(smarts_raw):
            smarts = smarts_raw
        else:
            st.sidebar.warning("Invalid SMARTS — substructure filter ignored.")

    # between() is inclusive on both ends, so a compound on a boundary (e.g. MW 500)
    # matches both the 450-500 and 500-550 windows.
    view = cat[
        cat["best_pchembl"].between(p_lo, p_hi)
        & cat["mw_freebase"].fillna(0).between(mw_lo, mw_hi)
        & cat["alogp"].fillna(0).between(logp_lo, logp_hi)
        & cat["hbd"].fillna(0).between(hbd_lo, hbd_hi)
        & cat["hba"].fillna(0).between(hba_lo, hba_hi)
        & (cat["num_ro5_violations"].fillna(0) <= max_ro5)
    ]
    if query:
        view = view[
            view["molecule_chembl_id"].str.lower().str.contains(query)
            | view["pref_name"].fillna("").str.lower().str.contains(query)
        ]
    if smarts:
        hits = _substructure_hits(smarts, tuple(view["canonical_smiles"].fillna("")))
        view = view[view["canonical_smiles"].fillna("").map(lambda s: hits.get(s, False))]
        st.caption(f"Substructure filter `{smarts}` active — matching atoms highlighted below.")

    disp = view.sort_values("best_pchembl", ascending=False).reset_index(drop=True)
    st.caption(f"{len(disp)} compounds — click a row to inspect it")
    # Tolerate a warehouse built before has_pdb existed: only show present columns.
    list_cols = [c for c in _LIST_COLS if c in disp.columns]
    # Mark the top compound on first open so the highlighted row matches the detail below.
    logic.preselect_first_row(st.session_state, "lib_rows")
    event = st.dataframe(
        disp[list_cols],
        hide_index=True,
        width="stretch",
        column_config=_list_column_config(),
        on_select="rerun",
        selection_mode="single-row",
        key="lib_rows",
    )
    if disp.empty:
        return

    # Selection is driven purely by clicking a row; default to the top compound.
    selected_rows = event.selection.rows
    idx = selected_rows[0] if selected_rows and selected_rows[0] < len(disp) else 0
    chosen = disp.iloc[idx]["molecule_chembl_id"]

    st.divider()
    fingerprints = data.load_fingerprints(con)
    compound_detail.render(
        con,
        disp.iloc[idx],
        chosen,
        fingerprints=fingerprints,
        catalog=catalog,
        highlight_smarts=smarts,
    )
