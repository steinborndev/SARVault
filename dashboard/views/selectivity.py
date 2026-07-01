"""Selectivity page (multi-target compounds only)."""

import streamlit as st

from dashboard import charts, data, logic


def render(con, scope):
    st.header("🎯 Selectivity")
    keys = logic.scope_compound_keys(data.load_target_sar(con), scope)
    sel = data.load_selectivity(con)
    sel = sel[sel["compound_key"].isin(keys)]
    multi = sel[sel["n_targets"] >= 2]
    st.caption(f"{len(multi)} multi-target compounds in scope (selectivity is defined only for these)")
    if multi.empty:
        st.info("No multi-target compounds in the current scope.")
        return
    st.plotly_chart(charts.selectivity_scatter(multi), use_container_width=True)
    st.dataframe(
        multi.sort_values("selectivity_index", ascending=False),
        hide_index=True,
        use_container_width=True,
    )
