"""Activity Cliffs page — the SAR centrepiece.

An activity cliff is a pair of structurally similar compounds whose potency
differs sharply *on the same target*: a small structural change with a large
biological consequence. This page lets the user tune the similarity and
Δ-potency thresholds, ranks pairs by SALI, and shows the two structures side by
side so the responsible change is visible.
"""

import base64

import pandas as pd
import streamlit as st

from dashboard import charts, chem, data


def _structure_img(smiles, highlight=None) -> str | None:
    svg = chem.smiles_to_svg(smiles, width=360, height=300, highlight_smarts=highlight)
    if not svg:
        return None
    b64 = base64.b64encode(svg.encode()).decode()
    return (
        f'<img src="data:image/svg+xml;base64,{b64}" '
        'style="width:100%; max-width:360px; background:#ffffff; border-radius:8px;">'
    )


def _cliff_column_config():
    return {
        "pair": st.column_config.TextColumn("Pair"),
        "target_pref_name": st.column_config.TextColumn("Target"),
        "tanimoto": st.column_config.ProgressColumn(
            "Tanimoto", format="%.3f", min_value=0.0, max_value=1.0
        ),
        "delta_pchembl": st.column_config.NumberColumn("Δ pChEMBL", format="%.2f"),
        "sali": st.column_config.NumberColumn("SALI", format="%.1f", help="|Δ| / (1 − Tanimoto)"),
        "same_scaffold": st.column_config.CheckboxColumn("Same scaffold"),
        "is_identical_fp": st.column_config.CheckboxColumn(
            "Identical 2D", help="Identical ECFP4 — likely stereo/tautomer/replicate, not a 2D change"
        ),
    }


def _pair_detail(row, smiles_by_id):
    a, b = row["molecule_chembl_id_a"], row["molecule_chembl_id_b"]
    st.markdown(f"#### {a}  ⇄  {b}  ·  {row['target_pref_name']}")

    m = st.columns(3)
    m[0].metric(f"{a} pChEMBL", f"{row['pchembl_a']:.2f}")
    m[1].metric(f"{b} pChEMBL", f"{row['pchembl_b']:.2f}")
    if pd.notna(row["sali"]):
        m[2].metric("SALI", f"{row['sali']:.1f}", f"Δ {row['delta_pchembl']:.2f} · T {row['tanimoto']:.3f}")
    else:
        m[2].metric("Δ pChEMBL", f"{row['delta_pchembl']:.2f}", "identical 2D fingerprint")

    left, right = st.columns(2)
    for col, cid in ((left, a), (right, b)):
        with col:
            img = _structure_img(smiles_by_id.get(cid))
            if img:
                st.markdown(img, unsafe_allow_html=True)
            else:
                st.info(f"No structure for {cid}.")
            st.caption(cid)

    if bool(row["is_identical_fp"]):
        st.caption(
            "These two share an identical 2D (ECFP4) fingerprint, so the potency "
            "gap comes from something the 2D graph doesn't capture — stereochemistry, "
            "tautomer/salt form, or measurement variance — rather than a structural edit."
        )


def render(con, scope):
    st.header("Activity cliffs")
    st.write(
        "Structurally similar compounds with a large potency difference on the same "
        "target — the sharpest signal in SAR. SALI = |Δ pChEMBL| / (1 − Tanimoto)."
    )

    cliffs = data.load_activity_cliffs(con)
    if cliffs.empty:
        st.info(
            "No activity-cliff mart in this warehouse (rebuild with the cheminfo + "
            "cliff models), or no cliffs in scope."
        )
        return

    targets = (scope or {}).get("targets")
    if targets:
        cliffs = cliffs[cliffs["target_pref_name"].isin(targets)]

    with st.sidebar:
        st.markdown("### Cliff thresholds")
        t_lo = float(cliffs["tanimoto"].min())
        min_tan = st.slider("Min Tanimoto", max(0.75, round(t_lo, 2)), 1.0, 0.80, 0.01)
        d_hi = float(max(4.0, cliffs["delta_pchembl"].max()))
        min_delta = st.slider("Min |Δ pChEMBL|", 1.0, round(d_hi, 1), 2.0, 0.5)
        show_identical = st.toggle(
            "Include identical-2D pairs", value=False,
            help="Tanimoto = 1 pairs (stereo/tautomer/replicate)",
        )

    from dashboard import logic

    view = logic.filter_cliffs(cliffs, min_tan, min_delta, show_identical)

    st.caption(
        f"{len(view)} cliff pairs at Tanimoto ≥ {min_tan:.2f} and |Δ pChEMBL| ≥ {min_delta:.1f}"
    )
    if view.empty:
        st.info("No cliffs at these thresholds — loosen the sliders in the sidebar.")
        return

    # SALI-ranked table drives both the scatter and the detail; _row ties them together.
    ranked = view.sort_values("sali", ascending=False, na_position="last").reset_index(drop=True)
    ranked["_row"] = range(len(ranked))
    ranked["pair"] = ranked["molecule_chembl_id_a"] + " ⇄ " + ranked["molecule_chembl_id_b"]

    scatter_event = st.plotly_chart(
        charts.cliff_scatter(ranked),
        width="stretch",
        on_select="rerun",
        selection_mode="points",
        key="cliff_scatter",
    )

    st.divider()
    list_cols = [
        "pair", "target_pref_name", "tanimoto", "delta_pchembl", "sali",
        "same_scaffold", "is_identical_fp",
    ]
    st.caption("Ranked by SALI — click a row, or click a point in the plot above.")
    # Mark the top-SALI pair on first open so the highlighted row matches the detail below.
    logic.preselect_first_row(st.session_state, "cliff_rows")
    table_event = st.dataframe(
        ranked[list_cols],
        hide_index=True,
        width="stretch",
        column_config=_cliff_column_config(),
        on_select="rerun",
        selection_mode="single-row",
        key="cliff_rows",
    )

    # Reconcile the two selectors: a click on either the scatter or the table wins over
    # a stale selection, and the last chosen row persists across unrelated reruns.
    scatter_idx = None
    if scatter_event.selection and scatter_event.selection.get("points"):
        cd = scatter_event.selection["points"][0].get("customdata")
        if cd:
            scatter_idx = int(cd[0])
    table_idx = table_event.selection.rows[0] if table_event.selection.rows else None

    prev = st.session_state.get("_cliff_sel", {"scatter": None, "table": None, "idx": 0})
    idx = prev["idx"]
    if scatter_idx is not None and scatter_idx != prev["scatter"]:
        idx = scatter_idx
    if table_idx is not None and table_idx != prev["table"]:
        idx = table_idx
    if idx is None or idx >= len(ranked):
        idx = 0
    st.session_state["_cliff_sel"] = {"scatter": scatter_idx, "table": table_idx, "idx": idx}

    # SMILES lookup for the side-by-side render (from the compound catalog).
    catalog = data.load_compound_catalog(con)
    smiles_by_id = dict(zip(catalog["molecule_chembl_id"], catalog["canonical_smiles"]))

    st.divider()
    _pair_detail(ranked.iloc[idx], smiles_by_id)
