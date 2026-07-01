"""Chemical-space page."""

import streamlit as st

from dashboard import charts, data, logic

_AXES = ["mw_freebase", "alogp", "psa", "qed_weighted", "rotatable_bonds"]


def render(con, scope):
    st.header("🧪 Chemical space")
    keys = logic.scope_compound_keys(data.load_target_sar(con), scope)
    chem = data.load_chemical_space(con)
    chem = chem[chem["compound_key"].isin(keys)]

    approved_only = st.checkbox("Approved drugs only", value=False)
    scoped = chem[chem["is_approved_drug"]] if approved_only else chem

    col_x, col_y = st.columns(2)
    x_col = col_x.selectbox("X axis", _AXES, index=0)
    y_col = col_y.selectbox("Y axis", _AXES, index=1)
    st.caption(f"{len(scoped)} compounds")
    st.plotly_chart(charts.chemical_space_scatter(scoped, x_col, y_col), use_container_width=True)
    st.plotly_chart(charts.property_histogram(scoped, x_col), use_container_width=True)
