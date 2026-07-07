"""Payload Classes: compare ADC-payload mechanism classes (tubulin vs topoisomerase)."""

import streamlit as st

from dashboard import charts, data, logic

_INTRO = """
Compounds grouped by **ADC-payload mechanism class** (from the target they act on),
so the tubulin and topoisomerase classes compare head to head. This is the chemistry
behind the ADC-payload shift from tubulin inhibitors (auristatins, maytansinoids)
toward topoisomerase-I inhibitors (camptothecin / exatecan).
"""


def render(con, scope):
    st.markdown(_INTRO)

    profile = data.load_payload_class_profile(con)
    if profile.empty:
        st.info("No payload-class profile is available in this warehouse yet.")
        return

    display = profile.copy()
    display["payload_class"] = logic.label_payload_class(display["payload_class"])
    display = display.rename(
        columns={
            "payload_class": "payload class",
            "n_compounds": "compounds",
            "n_measurements": "measurements",
            "median_pchembl": "median pChEMBL",
            "max_pchembl": "max pChEMBL",
            "p25_pchembl": "p25",
            "p75_pchembl": "p75",
            "n_sub_nanomolar": "sub-nM compounds",
        }
    )
    st.dataframe(display, hide_index=True, width="stretch")
    st.caption(
        "Per-class target potency (enzyme / binding pChEMBL). The potency distribution "
        "behind these numbers is on the SAR Ranking page (group by payload class)."
    )

    st.divider()
    st.subheader("Cellular cytotoxicity of reference payloads")
    cytotox = data.load_compound_cytotoxicity(con)
    if cytotox.empty:
        st.info("No cytotoxicity data is available in this warehouse yet.")
        return

    by_payload = logic.cytotox_by_payload(cytotox)
    chart_df = by_payload.assign(
        payload_class=logic.label_payload_class(by_payload["payload_class"])
    )
    st.plotly_chart(charts.cytotox_bar(chart_df), width="stretch")

    table = by_payload.rename(
        columns={
            "reference_name": "payload",
            "payload_class": "class",
            "best_p_cyto": "best cellular pChEMBL",
            "n_cell_lines": "cell lines",
        }
    )
    table["class"] = logic.label_payload_class(table["class"])
    st.dataframe(table, hide_index=True, width="stretch")
    st.caption(
        "Best cellular potency (p_cyto = -log10 of the GI50 / IC50 molar concentration) "
        "across tested cell lines, for the named clinical / reference ADC payloads. "
        "pChEMBL >= 9 is sub-nanomolar, the band that makes a viable payload."
    )
