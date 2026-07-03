"""Chemical Series page — compounds grouped into Bemis-Murcko scaffold series.

SAR is reasoned about per series, not per isolated compound. This page lists the
scaffolds present in the set with their size, potency spread and target reach, and
drills into a chosen series to show the shared scaffold and its member compounds.
"""

import base64

import pandas as pd
import streamlit as st

from dashboard import chem, compound_detail, data


def _structure_img(smiles, width: int = 320, height: int = 260) -> str | None:
    svg = chem.smiles_to_svg(smiles, width=width, height=height)
    if not svg:
        return None
    b64 = base64.b64encode(svg.encode()).decode()
    return (
        f'<img src="data:image/svg+xml;base64,{b64}" '
        'style="width:100%; max-width:320px; background:#ffffff; border-radius:8px;">'
    )


def _series_column_config():
    return {
        "murcko_scaffold_smiles": st.column_config.TextColumn("Scaffold (SMILES)"),
        "n_compounds": st.column_config.NumberColumn("Compounds", format="%d"),
        "n_targets": st.column_config.NumberColumn("Targets", format="%d"),
        "median_pchembl": st.column_config.NumberColumn("Median pChEMBL", format="%.2f"),
        "max_pchembl": st.column_config.NumberColumn("Max pChEMBL", format="%.2f"),
        "pchembl_range": st.column_config.NumberColumn(
            "Potency range", format="%.2f", help="Max − min best-pChEMBL across members"
        ),
        "top_compound": st.column_config.TextColumn("Most potent"),
    }


def _member_column_config():
    return {
        "molecule_chembl_id": st.column_config.TextColumn("ChEMBL ID"),
        "pref_name": st.column_config.TextColumn("Name"),
        "best_pchembl": st.column_config.ProgressColumn(
            "Best pChEMBL", format="%.2f", min_value=0.0, max_value=14.0
        ),
        "n_targets": st.column_config.NumberColumn("Targets", format="%d"),
    }


_CHEMBL_URL = "https://www.ebi.ac.uk/chembl/compound_report_card/{}/"


def render(con, scope):
    st.header("Chemical series")
    st.write(
        "Compounds grouped by their Bemis-Murcko scaffold — the chemical series a "
        "medicinal chemist reasons about. Each row is a scaffold shared by two or more "
        "compounds; a wide potency range within one series is where the SAR lives."
    )

    series = data.load_chemical_series(con)
    if series.empty:
        st.info(
            "No chemical-series mart in this warehouse (rebuild with the cheminfo + "
            "series models)."
        )
        return

    with st.sidebar:
        st.markdown("### Filters")
        max_size = int(series["n_compounds"].max())
        min_size = st.slider("Min compounds in series", 2, max(max_size, 2), 2)
        sort_by = st.selectbox(
            "Sort by",
            ["n_compounds", "pchembl_range", "max_pchembl", "n_targets"],
            format_func={
                "n_compounds": "series size",
                "pchembl_range": "potency range",
                "max_pchembl": "max potency",
                "n_targets": "target reach",
            }.get,
        )

    view = series[series["n_compounds"] >= min_size].sort_values(
        sort_by, ascending=False
    ).reset_index(drop=True)
    st.caption(f"{len(view)} multi-compound series — click a row to open it")

    list_cols = [
        "murcko_scaffold_smiles", "n_compounds", "n_targets",
        "median_pchembl", "max_pchembl", "pchembl_range", "top_compound",
    ]
    event = st.dataframe(
        view[list_cols],
        hide_index=True,
        width="stretch",
        column_config=_series_column_config(),
        on_select="rerun",
        selection_mode="single-row",
        key="series_rows",
    )
    if view.empty:
        return

    sel = event.selection.rows
    idx = sel[0] if sel and sel[0] < len(view) else 0
    row = view.iloc[idx]

    st.divider()
    left, right = st.columns([5, 7])
    with left:
        img = _structure_img(row["murcko_scaffold_smiles"])
        if img:
            st.markdown(img, unsafe_allow_html=True)
        st.caption(f"Scaffold: `{row['murcko_scaffold_smiles']}`")
    with right:
        m = st.columns(3)
        m[0].metric("Compounds", int(row["n_compounds"]))
        m[1].metric("Targets", int(row["n_targets"]))
        rng = row["pchembl_range"]
        m[2].metric("Potency range", f"{rng:.2f}" if pd.notna(rng) else "—")
        st.caption(
            f"Median best-pChEMBL {row['median_pchembl']:.2f} · most potent member "
            f"[{row['top_compound']}]({_CHEMBL_URL.format(row['top_compound'])})"
        )

    st.markdown("**Member compounds** — click a row to inspect it")
    members = data.scaffold_members(con, int(row["scaffold_key"]))
    if members.empty:
        return
    members_disp = members.drop(columns=["canonical_smiles"], errors="ignore")
    mevent = st.dataframe(
        members_disp,
        hide_index=True,
        width="stretch",
        column_config=_member_column_config(),
        on_select="rerun",
        selection_mode="single-row",
        key=f"members_{int(row['scaffold_key'])}",
    )
    msel = mevent.selection.rows
    midx = msel[0] if msel and msel[0] < len(members) else 0
    member = members.iloc[midx]

    st.divider()
    _member_context(member, midx, members, row)

    detail_row = data.compound_row(con, member["molecule_chembl_id"])
    if detail_row is None:
        st.info("No detailed record for this compound.")
        return
    compound_detail.render(
        con,
        detail_row,
        member["molecule_chembl_id"],
        fingerprints=data.load_fingerprints(con),
        catalog=data.load_compound_catalog(con),
    )


def _member_context(member, idx, members, series_row):
    """Where this compound sits within its scaffold series (rank + Δ to series stats)."""
    best = member["best_pchembl"]
    measured = int(members["best_pchembl"].notna().sum())
    cols = st.columns(3)
    if pd.notna(best):
        # members are sorted by best_pchembl desc (nulls last), so idx+1 is the rank.
        cols[0].metric("Potency rank in series", f"#{idx + 1} of {measured}")
        cols[1].metric(
            "Δ to series best", f"{best - series_row['max_pchembl']:+.2f}",
            help="This compound's best pChEMBL minus the series' most potent",
        )
        cols[2].metric(
            "Δ to series median", f"{best - series_row['median_pchembl']:+.2f}"
        )
    else:
        cols[0].metric("Potency rank in series", "—", "no activity data")
