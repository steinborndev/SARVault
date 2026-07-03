"""Chemical Series page - compounds grouped into Bemis-Murcko scaffold series.

SAR is reasoned about per series, not per isolated compound. This page lists the
scaffolds present in the set with their size, potency spread and target reach, and
drills into a chosen series to show the shared scaffold and its member compounds.
"""

import base64

import pandas as pd
import streamlit as st

from dashboard import chem, compound_detail, data, logic

# Enlarge the clickable rank count to the metric-value size so it matches the Delta
# metrics beside it. Scoped to the keyed nav container via Streamlit's documented
# `st-key-<key>` class, so it targets only these buttons and nothing else.
_RANK_NAV_CSS = """
<style>
[class*="st-key-sv_rank_nav"] button,
[class*="st-key-sv_rank_nav"] button * {
    font-size: 2.25rem !important;
    font-weight: 600 !important;
    line-height: 1.1 !important;
}
[class*="st-key-sv_rank_nav"] button { padding: 0 0.35rem; min-height: 0; }
[class*="st-key-sv_rank_nav"] [data-testid="stMarkdownContainer"] p {
    font-size: 2.25rem;
    line-height: 1.1;
    margin: 0;
}
</style>
"""


@st.cache_data(show_spinner=False)
def _series_frame(scaffold_smiles, member_smiles):
    """Cached shared drawing window for a series (recomputed only when it changes)."""
    return chem.scaffold_frame(scaffold_smiles, list(member_smiles))


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
        "Compounds grouped by their Bemis-Murcko scaffold - the chemical series a "
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
    st.caption(f"{len(view)} multi-compound series - click a row to open it")

    list_cols = [
        "murcko_scaffold_smiles", "n_compounds", "n_targets",
        "median_pchembl", "max_pchembl", "pchembl_range", "top_compound",
    ]
    # Mark the top series on first open so the highlighted row matches the detail below.
    logic.preselect_first_row(st.session_state, "series_rows")
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
        m[2].metric("Potency range", f"{rng:.2f}" if pd.notna(rng) else "-")
        st.caption(
            f"Median best-pChEMBL {row['median_pchembl']:.2f} · most potent member "
            f"[{row['top_compound']}]({_CHEMBL_URL.format(row['top_compound'])})"
        )

    st.divider()
    st.markdown("**Member compounds** - click a row to inspect it")
    members = data.scaffold_members(con, int(row["scaffold_key"]))
    if members.empty:
        return
    members_disp = members.drop(columns=["canonical_smiles"], errors="ignore")
    # The member table is keyed per scaffold, so opening a new series produces a fresh
    # key - seeding row 0 marks the first member of whichever series is being viewed.
    members_key = f"members_{int(row['scaffold_key'])}"
    logic.preselect_first_row(st.session_state, members_key)
    mevent = st.dataframe(
        members_disp,
        hide_index=True,
        width="stretch",
        column_config=_member_column_config(),
        on_select="rerun",
        selection_mode="single-row",
        key=members_key,
    )
    msel = mevent.selection.rows
    midx = msel[0] if msel and msel[0] < len(members) else 0
    member = members.iloc[midx]

    st.divider()
    _member_context(member, midx, members, row, members_key)

    dcol = st.columns([1, 1, 4])
    align = dcol[0].toggle(
        "Align to scaffold", value=True, key="series_align",
        help="Draw every member in the scaffold's orientation, so the shared core stays "
        "fixed and only the substituents move as you step through the series.",
    )
    highlight = dcol[1].toggle(
        "Highlight core", value=False, key="series_highlight",
        help="Wash the shared scaffold with a subtle tint so substituent changes stand out.",
    )

    detail_row = data.compound_row(con, member["molecule_chembl_id"])
    if detail_row is None:
        st.info("No detailed record for this compound.")
        return
    # Shared drawing window so the aligned core is pixel-stable across the whole series
    # (computed once per series and cached, not per step). Only needed when aligning.
    frame = (
        _series_frame(row["murcko_scaffold_smiles"], tuple(members["canonical_smiles"].fillna("")))
        if align
        else None
    )
    compound_detail.render(
        con,
        detail_row,
        member["molecule_chembl_id"],
        fingerprints=data.load_fingerprints(con),
        catalog=data.load_compound_catalog(con),
        scaffold_smiles=row["murcko_scaffold_smiles"],
        align_to_scaffold=align,
        highlight_scaffold=highlight,
        frame=frame,
    )


def _member_context(member, idx, members, series_row, members_key):
    """Rank within the scaffold series plus Delta to series stats.

    The rank count doubles as the navigator (no separate arrow controls): clicking the
    position number steps to the previous (more potent) member, clicking the total steps
    to the next, by moving the dataframe's selection via ``logic.step_selection``.
    """
    best = member["best_pchembl"]
    n = len(members)
    cols = st.columns([3, 4, 4], vertical_alignment="bottom")
    with cols[0]:
        st.caption("Potency rank in series")
        st.markdown(_RANK_NAV_CSS, unsafe_allow_html=True)
        with st.container(
            horizontal=True, gap="small", vertical_alignment="center", key="sv_rank_nav"
        ):
            st.button(
                f"#{idx + 1}", key=f"prev_{members_key}", type="tertiary",
                disabled=idx <= 0, help="Previous (more potent) member",
                on_click=logic.step_selection, args=(st.session_state, members_key, -1, n),
            )
            st.markdown("of")
            st.button(
                f"{n}", key=f"next_{members_key}", type="tertiary",
                disabled=idx >= n - 1, help="Next (less potent) member",
                on_click=logic.step_selection, args=(st.session_state, members_key, 1, n),
            )
    if pd.notna(best):
        cols[1].metric(
            "Δ to series best", f"{best - series_row['max_pchembl']:+.2f}",
            help="This compound's best pChEMBL minus the series' most potent",
        )
        cols[2].metric("Δ to series median", f"{best - series_row['median_pchembl']:+.2f}")
