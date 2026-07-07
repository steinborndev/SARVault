"""Landing page: project intro, KPIs, per-target overview and page guide."""

import streamlit as st

from dashboard import data, logic

_GUIDE = """
**Use the pages in the sidebar:**

- **Compound Library** — browse and filter all compounds; open one for structure, properties, its per-target potency and structural analogs.
- **SAR Ranking** — rank compounds by median potency for a chosen target.
- **Activity Cliffs** — pairs of similar compounds with a large potency gap on the same target, ranked by SALI.
- **Chemical Series** — compounds grouped by Bemis-Murcko scaffold; series size, potency spread and target reach.
- **Selectivity** — multi-target compounds: selectivity index (best vs. second-best target) against potency.
- **Chemical Space** — physicochemical profile with an approved-vs-research view.
- **Data Quality** — provenance: ChEMBL release, rows per layer, assay-confidence distribution.

The **Scope** selector (top right) narrows every page by target, approval status and minimum potency.
"""


def render(con, scope):
    st.header("Overview")
    st.caption("The warehouse at a glance: scale, coverage and a guide to the pages.")
    target_sar = data.load_target_sar(con)
    catalog = data.load_compound_catalog(con)
    metrics = logic.overview_metrics(target_sar, catalog, scope)

    row1 = st.columns(3)
    row1[0].metric("Compounds", f"{metrics['compounds']:,}")
    row1[1].metric("Activities (fact)", f"{metrics['activities']:,}")
    row1[2].metric("Targets", metrics["targets"])
    row2 = st.columns(3)
    row2[0].metric("Compound-target pairs", f"{metrics['pairs']:,}")
    row2[1].metric("Multi-target compounds", metrics["multi_target"])
    row2[2].metric("Approved drugs", metrics["approved"])

    st.divider()
    st.subheader("Per-target overview")
    summary = data.target_summary(con)
    targets = (scope or {}).get("targets")
    if targets:
        summary = summary[summary["target"].isin(targets)]
    st.dataframe(summary, hide_index=True, width="stretch")

    st.divider()
    st.markdown(_GUIDE)
    st.caption(
        "Source: ChEMBL (EMBL-EBI), released under CC BY-SA. Read-only view over the "
        "dbt-modelled warehouse - see LICENSE-DATA.md."
    )
