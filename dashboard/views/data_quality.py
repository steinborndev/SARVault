"""Data quality & provenance page (warehouse-wide)."""

import streamlit as st

from dashboard import charts, data


def render(con, scope):  # scope intentionally unused: this view is warehouse-wide
    st.header("Data quality & provenance")

    cfg = data.pipeline_config()
    cols = st.columns(3)
    cols[0].metric("ChEMBL release", cfg["chembl_version"])
    cols[1].metric("Confidence floor", cfg["min_confidence_score"])
    cols[2].metric("Organism", cfg["organism"])

    st.divider()
    st.subheader("Rows per layer")
    st.dataframe(data.layer_counts(con), hide_index=True, width="stretch")

    st.divider()
    st.subheader("Assay confidence distribution")
    conf = data.confidence_distribution(con)
    st.plotly_chart(charts.confidence_bar(conf), width="stretch")
    st.dataframe(conf, hide_index=True, width="stretch")

    st.divider()
    st.subheader("Endpoint (standard type) distribution")
    stype = data.standard_type_distribution(con)
    st.plotly_chart(charts.standard_type_bar(stype), width="stretch")
    st.dataframe(stype, hide_index=True, width="stretch")

    st.divider()
    st.subheader("Per-target breakdown")
    st.dataframe(data.target_summary(con), hide_index=True, width="stretch")

    st.caption("Source: ChEMBL (EMBL-EBI), CC BY-SA. See LICENSE-DATA.md.")
