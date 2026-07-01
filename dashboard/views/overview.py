"""Landing page: project intro, KPIs, per-target overview and page guide."""

import streamlit as st

from dashboard import data, logic

_GUIDE = """
**Use the pages in the sidebar:**

- **📊 SAR Ranking** — rank compounds by median potency for a chosen target, with measurement and pChEMBL filters.
- **🎯 Selectivity** — multi-target compounds only: selectivity index (best vs. second-best target) against potency.
- **🧪 Chemical Space** — physicochemical profile with an approved-vs-research view over the payload chemical space.
- **🔎 Data Quality** — provenance: ChEMBL release, rows per layer, and the assay-confidence distribution.

The **🔬 Scope** selector in the sidebar narrows every page to the selected targets.
"""


def render(con, scope):
    st.subheader("Payload SAR Warehouse")
    st.write(
        "A layered warehouse over public **ChEMBL** bioactivity data, scoped to "
        "cytotoxic / tubulin-targeting compounds — the chemistry behind ADC payloads "
        "and classical chemotherapeutics."
    )

    sar = data.load_target_sar(con)
    metrics = logic.overview_metrics(
        sar, data.load_selectivity(con), data.load_chemical_space(con), scope
    )

    row1 = st.columns(3)
    row1[0].metric("Compounds", f"{metrics['compounds']:,}")
    row1[1].metric("Activities (fact)", f"{metrics['activities']:,}")
    row1[2].metric("Targets", metrics["targets"])
    row2 = st.columns(3)
    row2[0].metric("Compound-target pairs", f"{metrics['pairs']:,}")
    row2[1].metric("Multi-target compounds", metrics["multi_target"])
    row2[2].metric("Approved drugs", metrics["approved"])

    st.subheader("Per-target overview")
    summary = data.target_summary(con)
    if scope:
        summary = summary[summary["target"].isin(scope)]
    st.dataframe(summary, hide_index=True, use_container_width=True)

    st.markdown(_GUIDE)
    st.caption(
        "Source: ChEMBL (EMBL-EBI), released under CC BY-SA. Read-only view over the "
        "dbt-modelled warehouse - see LICENSE-DATA.md."
    )
