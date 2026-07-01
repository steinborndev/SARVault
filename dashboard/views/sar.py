"""SAR ranking page."""

import streamlit as st

from dashboard import charts, data, logic


def render(con, scope):
    st.header("📊 SAR ranking")
    target_sar = data.load_target_sar(con)
    catalog = data.load_compound_catalog(con)
    keys = logic.resolve_scope_keys(target_sar, catalog, scope)
    sar = logic.scoped_target_sar(target_sar, scope, keys)
    if sar.empty:
        st.info("No data in the current scope.")
        return

    targets = sorted(sar["target_pref_name"].unique())
    selected = st.selectbox("Target", targets)
    max_meas = int(sar["n_measurements"].max())

    with st.sidebar:
        st.markdown("### Filters")
        min_meas = st.slider("Min measurements", 1, max_meas, 1) if max_meas > 1 else 1
        low, high = st.slider("median pChEMBL range", 0.0, 14.0, (5.0, 14.0), step=0.5)

    view = sar[
        (sar["target_pref_name"] == selected)
        & (sar["n_measurements"] >= min_meas)
        & (sar["median_pchembl"].between(low, high))
    ]
    st.caption(f"{len(view)} compound-target pairs")
    st.plotly_chart(charts.sar_ranking_bar(view), width="stretch")
    st.dataframe(
        view.sort_values("median_pchembl", ascending=False),
        hide_index=True,
        width="stretch",
    )
