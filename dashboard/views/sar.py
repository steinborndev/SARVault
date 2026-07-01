"""SAR ranking page."""

import streamlit as st

from dashboard import charts, data, logic


def render(con, scope):
    st.header("📊 SAR ranking")
    sar = logic.scoped_target_sar(data.load_target_sar(con), scope)
    if sar.empty:
        st.info("No data in the current scope.")
        return

    targets = sorted(sar["target_pref_name"].unique())
    selected = st.selectbox("Target", targets)
    max_meas = int(sar["n_measurements"].max())
    min_meas = st.slider("Min measurements", 1, max(max_meas, 1), 1)
    low, high = st.slider("median pChEMBL range", 0.0, 14.0, (5.0, 14.0), 0.1)

    view = sar[
        (sar["target_pref_name"] == selected)
        & (sar["n_measurements"] >= min_meas)
        & (sar["median_pchembl"].between(low, high))
    ]
    st.caption(f"{len(view)} compound-target pairs")
    st.plotly_chart(charts.sar_ranking_bar(view), use_container_width=True)
    st.dataframe(
        view.sort_values("median_pchembl", ascending=False),
        hide_index=True,
        use_container_width=True,
    )
