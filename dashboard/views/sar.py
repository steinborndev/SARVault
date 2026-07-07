"""SAR ranking page."""

import streamlit as st

from dashboard import charts, data, logic


def render(con, scope):
    st.header("SAR ranking")
    st.caption(
        "Rank compounds by potency, and see how potency is distributed across targets "
        "or payload classes."
    )
    target_sar = data.load_target_sar(con)
    catalog = data.load_compound_catalog(con)
    keys = logic.resolve_scope_keys(target_sar, catalog, scope)
    sar = logic.scoped_target_sar(target_sar, scope, keys)
    if sar.empty:
        st.info("No data in the current scope.")
        return

    st.subheader("Potency landscape")
    # The Payload class grouping needs the payload_class column (added to
    # mart_target_sar in a later build); only offer the toggle when it is present,
    # so an older warehouse shows the target view rather than a dead control.
    has_class = "payload_class" in sar.columns
    group = (
        st.radio("Group by", ["Target", "Payload class"], horizontal=True, key="sar_violin_group")
        if has_class
        else "Target"
    )
    if group == "Payload class":
        grouped = sar.assign(payload_class=logic.label_payload_class(sar["payload_class"]))
        st.plotly_chart(
            charts.target_potency_violin(
                grouped, group_col="payload_class", group_label="payload class"
            ),
            width="stretch",
        )
        st.caption(
            "Each violin is the distribution of per-compound-target median pChEMBL for one "
            "ADC-payload mechanism class."
        )
    else:
        st.plotly_chart(charts.target_potency_violin(sar), width="stretch")
        st.caption(
            "Each violin is the distribution of per-compound median pChEMBL for one target."
        )

    st.divider()
    st.subheader("Rank compounds for one target")
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
