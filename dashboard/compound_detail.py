"""Shared per-compound detail card (structure, physchem, potency, analogs, 3D, Ro5).

Rendered from both the Compound Library and the Chemical Series pages so the two
show an identical, rich drill-down. ``row`` is a compound record with the
dim_compound / mart_compound_catalog columns (see ``data.compound_row``).
"""

import base64

import pandas as pd
import streamlit as st
import streamlit.components.v1 as components

from dashboard import charts, chem, data, logic, viewer

CHEMBL_URL = "https://www.ebi.ac.uk/chembl/compound_report_card/{}/"

_SLIM_PROPS = [
    ("HBA", "hba"),
    ("HBD", "hbd"),
    ("Rotatable bonds", "rotatable_bonds"),
    ("Aromatic rings", "aromatic_rings"),
    ("Ro3 pass", "ro3_pass"),
    ("Max phase", "max_phase"),
]

_METHOD_SHORT = {
    "X-ray diffraction": "X-ray",
    "Electron Microscopy": "cryo-EM",
    "Solution NMR": "NMR",
    "Solid-state NMR": "ss-NMR",
    "Neutron Diffraction": "neutron",
}


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


def _pdb_meta(rec) -> str:
    method = rec.get("method")
    method_s = _METHOD_SHORT.get(method, method) if isinstance(method, str) and method else None
    res = rec.get("resolution")
    res_s = f"{float(res):.1f} Å" if pd.notna(res) else None
    year = rec.get("year")
    year_s = str(int(year)) if pd.notna(year) else None
    method_res = " ".join(part for part in [method_s, res_s] if part)
    return " · ".join(part for part in [method_res, year_s] if part)


def _pdb_label(rec) -> str:
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


def _analog_column_config():
    return {
        "molecule_chembl_id": st.column_config.LinkColumn(
            "ChEMBL ID",
            help="Open the compound report card on ChEMBL",
            display_text=r"(CHEMBL\d+)",
        ),
        "pref_name": st.column_config.TextColumn("Name"),
        "tanimoto": st.column_config.ProgressColumn(
            "Tanimoto", format="%.3f", min_value=0.0, max_value=1.0
        ),
        "best_pchembl": st.column_config.NumberColumn("Best pChEMBL", format="%.2f"),
        "delta_pchembl": st.column_config.NumberColumn(
            "Δ pChEMBL", format="%.2f", help="Analog best pChEMBL − this compound's"
        ),
        "best_target": st.column_config.TextColumn("Best target"),
    }


def render(con, row, chosen, fingerprints=None, catalog=None, highlight_smarts=None):
    """Render the full compound detail card for ``row`` (identified by ``chosen``)."""
    name = row.get("pref_name")
    has_name = isinstance(name, str) and name.strip() and name != chosen
    title = f"{chosen} — {name}" if has_name else chosen
    suffix = " :green[· approved]" if bool(row.get("is_approved_drug")) else ""

    st.markdown(f"### {title}{suffix}")

    left, right = st.columns([5, 7])
    with left:
        svg = chem.smiles_to_svg(row.get("canonical_smiles"), highlight_smarts=highlight_smarts)
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
        elif len(profile) > 1:
            st.plotly_chart(charts.compound_potency_bar(profile), width="stretch")
            st.dataframe(profile, hide_index=True, width="stretch")
            if pd.notna(row.get("selectivity_index")):
                st.metric(
                    "Selectivity index (log10 fold)",
                    round(float(row["selectivity_index"]), 2),
                )
        else:
            st.caption("No measured activity for this compound.")
    with ref_col:
        st.markdown("**External references**")
        parts = [f"[ChEMBL]({CHEMBL_URL.format(chosen)})"]
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
        with st.expander(f"3D co-crystal structure (PDBe) — {len(pdb)} entries"):
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

    # --- structural analogs (ECFP4 Tanimoto over the fingerprint mart) ---
    ckey = int(row["compound_key"])
    has_fp = (
        fingerprints is not None
        and catalog is not None
        and not fingerprints.empty
        and (fingerprints["compound_key"] == ckey).any()
    )
    if has_fp:
        st.divider()
        st.markdown("**Structural analogs** — nearest compounds by ECFP4 Tanimoto")
        ctrl = st.columns(2)
        top_n = ctrl[0].slider("How many", 3, 25, 10, key=f"nn_n_{chosen}")
        min_sim = ctrl[1].slider("Min Tanimoto", 0.0, 1.0, 0.30, 0.05, key=f"nn_s_{chosen}")
        analogs = logic.nearest_neighbors(
            ckey, fingerprints, catalog, top_n=top_n, min_similarity=min_sim
        )
        if analogs.empty:
            st.info("No analogs above this Tanimoto threshold in the warehouse.")
        else:
            linked = analogs.copy()
            linked["molecule_chembl_id"] = linked["molecule_chembl_id"].map(CHEMBL_URL.format)
            st.dataframe(
                linked, hide_index=True, width="stretch",
                column_config=_analog_column_config(),
            )
            st.caption(
                "Δ pChEMBL is the analog's best potency minus this compound's — positive "
                "means a more potent close neighbour (a lead for the SAR series)."
            )

    st.divider()
    st.markdown(f"**Lipinski Ro5** — {_ro5_line(row)}")
    st.caption(
        "Ro5 is only weakly predictive here: ADC payloads / cytotoxics are "
        "antibody-delivered rather than orally absorbed, and highly potent payloads "
        "often violate it by design."
    )
