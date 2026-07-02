"""Compound library: filterable list + per-compound detail drill-down."""

import base64
import math

import pandas as pd
import streamlit as st
import streamlit.components.v1 as components

from dashboard import charts, chem, data, logic, viewer

_CHEMBL_URL = "https://www.ebi.ac.uk/chembl/compound_report_card/{}/"

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
        "has_pdb": st.column_config.CheckboxColumn("🧬 3D", help="Has a co-crystal PDB entry"),
        "best_pchembl": st.column_config.ProgressColumn(
            "Best pChEMBL", format="%.2f", min_value=0.0, max_value=14.0
        ),
        "best_target": st.column_config.TextColumn("Best target"),
        "n_targets": st.column_config.NumberColumn("Targets", format="%d"),
        "mw_freebase": st.column_config.NumberColumn("MW", format="%.0f"),
        "alogp": st.column_config.NumberColumn("logP", format="%.1f"),
        "num_ro5_violations": st.column_config.NumberColumn("Ro5 viol.", format="%d"),
    }

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


_METHOD_SHORT = {
    "X-ray diffraction": "X-ray",
    "Electron Microscopy": "cryo-EM",
    "Solution NMR": "NMR",
    "Solid-state NMR": "ss-NMR",
    "Neutron Diffraction": "neutron",
}


def _pdb_meta(rec) -> str:
    """Compact 'method resolution · year' string from a PDB entry row (parts optional)."""
    method = rec.get("method")
    method_s = _METHOD_SHORT.get(method, method) if isinstance(method, str) and method else None
    res = rec.get("resolution")
    res_s = f"{float(res):.1f} Å" if pd.notna(res) else None
    year = rec.get("year")
    year_s = str(int(year)) if pd.notna(year) else None
    method_res = " ".join(part for part in [method_s, res_s] if part)
    return " · ".join(part for part in [method_res, year_s] if part)


def _pdb_label(rec) -> str:
    """Dropdown label: 'CODE — full title · X-ray 2.1 Å · 2013' (metadata parts optional)."""
    pid = str(rec["pdb_id"]).upper()
    title = rec.get("title")
    head = f"{pid} — {title}" if isinstance(title, str) and title else pid
    meta = _pdb_meta(rec)
    return f"{head} · {meta}" if meta else head


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

    st.markdown(f"### {title}{suffix}")

    # Structure-led split: molecule on the left, all physicochemical readouts
    # (key metrics + slim property table) grouped together on the right.
    left, right = st.columns([5, 7])
    with left:
        svg = chem.smiles_to_svg(row.get("canonical_smiles"))
        if svg:
            st.markdown(_structure_html(svg), unsafe_allow_html=True)
        else:
            st.info("No structure available.")
    with right:
        metrics = st.columns(4)
        metrics[0].metric("MW", _mfmt(row["mw_freebase"], 0))
        metrics[1].metric("logP", _mfmt(row["alogp"], 1))
        metrics[2].metric("TPSA", _mfmt(row["psa"], 0))
        metrics[3].metric("QED", _mfmt(row["qed_weighted"], 2))
        slim = pd.DataFrame(
            [{"property": label, "value": _fmt(row.get(col))} for label, col in _SLIM_PROPS]
        )
        st.dataframe(slim, hide_index=True, width="stretch")

    # Potency and external references sit side by side beneath the split.
    pot_col, ref_col = st.columns(2)
    with pot_col:
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
                st.metric(
                    "Selectivity index (log10 fold)",
                    round(float(row["selectivity_index"]), 2),
                )
    with ref_col:
        st.markdown("**External references**")
        # ChEMBL always leads the list (replaces the old "View on ChEMBL" button).
        parts = [f"[ChEMBL]({_CHEMBL_URL.format(chosen)})"]
        xrefs = data.compound_xrefs(con, int(row["compound_key"]))
        for _, ref in xrefs.iterrows():
            extra = f" :gray[+{int(ref['n_refs']) - 1}]" if ref["n_refs"] > 1 else ""
            if pd.notna(ref["url"]) and ref["url"]:
                parts.append(f"[{ref['display_name']}]({ref['url']}){extra}")
            else:
                parts.append(f"{ref['display_name']} ({ref['xref_id']}){extra}")
        st.markdown(" · ".join(parts))

    pdb = data.compound_pdb_entries(con, int(row["compound_key"]))
    if not pdb.empty:
        with st.expander(f"🧬 3D co-crystal structure (PDBe) — {len(pdb)} entries"):
            labels = {rec["pdb_id"]: _pdb_label(rec) for _, rec in pdb.iterrows()}
            pick = st.selectbox(
                "PDB entry",
                pdb["pdb_id"].tolist(),
                format_func=lambda p: labels.get(p, p.upper()),
                key=f"pdb_pick_{chosen}",
            )
            sel = pdb.loc[pdb["pdb_id"] == pick].iloc[0]
            st.caption(
                f"Ligand `{sel['ligand_code']}` · "
                f"[open {pick.upper()} on PDBe](https://www.ebi.ac.uk/pdbe/entry/pdb/{pick})"
            )
            components.html(viewer.pdbe_molstar_html(pick), height=960)

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
    # Tolerate a warehouse built before has_pdb existed: only show present columns.
    list_cols = [c for c in _LIST_COLS if c in disp.columns]
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

    selected_rows = event.selection.rows
    if selected_rows and selected_rows[0] < len(disp):
        st.session_state["inspect_compound"] = disp.iloc[selected_rows[0]]["molecule_chembl_id"]

    st.subheader("Inspect a compound")
    options = disp["molecule_chembl_id"].tolist()
    if st.session_state.get("inspect_compound") not in options:
        st.session_state["inspect_compound"] = options[0]
    chosen = st.selectbox("Compound", options, key="inspect_compound", label_visibility="collapsed")
    _detail(con, disp[disp["molecule_chembl_id"] == chosen].iloc[0], chosen)
