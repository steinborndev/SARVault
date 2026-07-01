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

_SLIM_PROPS = [
    ("HBA", "hba"),
    ("HBD", "hbd"),
    ("Rotatable bonds", "rotatable_bonds"),
    ("Aromatic rings", "aromatic_rings"),
    ("Ro3 pass", "ro3_pass"),
    ("Max phase", "max_phase"),
]


def _fmt(value) -> str:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return "—"
    if isinstance(value, float):
        return f"{value:.2f}"
    return str(value)


def _num(value) -> str:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return "—"
    fv = round(float(value), 1)
    return str(int(fv)) if fv == int(fv) else f"{fv:.1f}"


def _mfmt(value, decimals: int) -> str:
    return f"{value:.{decimals}f}" if pd.notna(value) else "—"


def _structure_html(svg: str) -> str:
    b64 = base64.b64encode(svg.encode()).decode()
    return (
        f'<img src="data:image/svg+xml;base64,{b64}" '
        'style="width:100%; max-width:440px; background:#ffffff; border-radius:8px;">'
    )


def _int_max(series, default: int) -> int:
    clean = series.dropna()
    return int(clean.max()) if not clean.empty else default


def _ro5_line(row) -> str:
    breakdown = logic.ro5_breakdown(row)
    marks = []
    for item in breakdown["items"]:
        mark = "✓" if item["pass"] else ("✗" if item["pass"] is False else "—")
        marks.append(f"{item['label']} {mark} ({_num(item['value'])})")
    line = " · ".join(marks) + f" → {breakdown['violations']} violation(s)"
    official = row.get("num_ro5_violations")
    if pd.notna(official):
        line += f" · ChEMBL: {int(official)}"
    return line


def _detail(con, row, chosen):
    name = row.get("pref_name")
    title = f"{chosen} — {name}" if name and str(name) != chosen else chosen
    suffix = " :green[· approved]" if bool(row["is_approved_drug"]) else ""

    head_l, head_r = st.columns([3, 1])
    head_l.markdown(f"### {title}{suffix}")
    head_r.link_button("View on ChEMBL", _CHEMBL_URL.format(chosen))

    left, right = st.columns([3, 2])
    with left:
        svg = chem.smiles_to_svg(row.get("canonical_smiles"))
        if svg:
            st.markdown(_structure_html(svg), unsafe_allow_html=True)
        else:
            st.info("No structure available.")
    with right:
        top = st.columns(2)
        top[0].metric("MW", _mfmt(row["mw_freebase"], 0))
        top[1].metric("logP", _mfmt(row["alogp"], 1))
        bottom = st.columns(2)
        bottom[0].metric("TPSA", _mfmt(row["psa"], 0))
        bottom[1].metric("QED", _mfmt(row["qed_weighted"], 2))
        slim = pd.DataFrame(
            [{"property": label, "value": _fmt(row.get(col))} for label, col in _SLIM_PROPS]
        )
        st.dataframe(slim, hide_index=True, width="stretch")

    st.markdown("**Per-target potency**")
    profile = data.compound_target_profile(con, int(row["compound_key"]))
    if len(profile) == 1:
        p = profile.iloc[0]
        st.markdown(
            f"**{p['target']}** — median pChEMBL {p['median_pchembl']:.2f} · "
            f"max {p['max_pchembl']:.2f} · {int(p['n_measurements'])} measurement(s)"
        )
    else:
        st.plotly_chart(charts.compound_potency_bar(profile), width="stretch")
        st.dataframe(profile, hide_index=True, width="stretch")
        if pd.notna(row["selectivity_index"]):
            st.metric("Selectivity index (log10 fold)", round(float(row["selectivity_index"]), 2))

    st.divider()
    st.markdown(f"**Lipinski Ro5** — {_ro5_line(row)}")
    st.caption(
        "Ro5 is only weakly predictive here: ADC payloads / cytotoxics are "
        "antibody-delivered rather than orally absorbed, and highly potent payloads "
        "often violate it by design."
    )


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
        hbd_max = _int_max(cat["hbd"], 10)
        hbd_lo, hbd_hi = st.slider("HBD range", 0, max(hbd_max, 1), (0, max(hbd_max, 1)))
        hba_max = _int_max(cat["hba"], 15)
        hba_lo, hba_hi = st.slider("HBA range", 0, max(hba_max, 1), (0, max(hba_max, 1)))
        max_ro5 = st.slider("Max Ro5 violations", 0, 4, 4)

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

    selected_rows = event.selection.rows
    if selected_rows and selected_rows[0] < len(disp):
        st.session_state["inspect_compound"] = disp.iloc[selected_rows[0]]["molecule_chembl_id"]

    st.subheader("Inspect a compound")
    options = disp["molecule_chembl_id"].tolist()
    if st.session_state.get("inspect_compound") not in options:
        st.session_state["inspect_compound"] = options[0]
    chosen = st.selectbox("Compound", options, key="inspect_compound", label_visibility="collapsed")
    _detail(con, disp[disp["molecule_chembl_id"] == chosen].iloc[0], chosen)
