"""Selectivity page (multi-target compounds only)."""

import streamlit as st

from dashboard import charts, data, logic


def render(con, scope):
    st.header("🎯 Selectivity")
    target_sar = data.load_target_sar(con)
    catalog = data.load_compound_catalog(con)
    keys = logic.resolve_scope_keys(target_sar, catalog, scope)
    sel = data.load_selectivity(con)
    sel = sel[sel["compound_key"].isin(keys)]
    multi = sel[sel["n_targets"] >= 2]
    st.caption(f"{len(multi)} multi-target compounds in scope (selectivity is defined only for these)")
    if multi.empty:
        st.info("No multi-target compounds in the current scope.")
        return
    st.plotly_chart(charts.selectivity_scatter(multi), width="stretch")
    st.dataframe(
        multi.sort_values("selectivity_index", ascending=False),
        hide_index=True,
        width="stretch",
    )

    st.subheader("Selectivity matrix")
    sar = logic.scoped_target_sar(target_sar, scope, keys)
    heat = sar[sar["compound_key"].isin(multi["compound_key"])]
    st.plotly_chart(charts.sar_heatmap(heat), width="stretch")
    st.caption(
        "Median pChEMBL per target for the most potent multi-target compounds — "
        "a single bright cell reads as selective, several as promiscuous."
    )
