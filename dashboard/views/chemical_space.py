"""Chemical-space page."""

import streamlit as st

from dashboard import charts, data, logic
from dashboard.chem import heavy_atom_count

_AXES = ["mw_freebase", "alogp", "psa", "qed_weighted", "rotatable_bonds"]


@st.cache_data(show_spinner=False)
def _hac_by_smiles(smiles: tuple[str, ...]) -> dict:
    """Heavy-atom count per unique SMILES (cached; RDKit parse is the cost)."""
    return {s: heavy_atom_count(s) for s in smiles}


def render(con, scope):
    st.header("Chemical space")
    st.caption(
        "The physicochemical profile of the compound set, with an "
        "approved-vs-research lens."
    )
    target_sar = data.load_target_sar(con)
    catalog = data.load_compound_catalog(con)
    keys = logic.resolve_scope_keys(target_sar, catalog, scope)
    chem = data.load_chemical_space(con)
    chem = chem[chem["compound_key"].isin(keys)]
    if chem.empty:
        st.info("No compounds in the current scope.")
        return

    st.caption(f"{len(chem)} compounds")

    has_embedding = "umap_x" in chem.columns and chem["umap_x"].notna().any()
    modes = ["Structural embedding (UMAP)", "Property axes"] if has_embedding else ["Property axes"]
    mode = st.radio("View", modes, horizontal=True)

    if mode == "Structural embedding (UMAP)":
        st.caption(
            "2-D UMAP of the ECFP4 fingerprints — proximity ≈ structural similarity. "
            "Colour reveals where potency (or an approval class, or a scaffold series) "
            "concentrates in chemical space."
        )
        color_by = st.selectbox(
            "Colour by", ["potency", "approval", "series"],
            format_func={"potency": "best potency", "approval": "approval status",
                         "series": "scaffold series"}.get,
        )
        st.plotly_chart(charts.embedding_scatter(chem, color_by), width="stretch")
    else:
        col_x, col_y = st.columns(2)
        x_col = col_x.selectbox("X axis", _AXES, index=0)
        y_col = col_y.selectbox("Y axis", _AXES, index=1)
        st.plotly_chart(charts.chemical_space_scatter(chem, x_col, y_col), width="stretch")
        st.divider()
        st.plotly_chart(charts.property_histogram(chem, x_col), width="stretch")

    st.divider()
    st.subheader("Ligand efficiency")
    st.caption(
        "LE = 1.37 × pChEMBL / heavy-atom count (potency per atom); "
        "LLE = pChEMBL − logP (potency vs. lipophilicity). For ADC payloads, "
        "high efficiency at high potency (upper right) is the sweet spot."
    )
    smiles = catalog[["compound_key", "canonical_smiles"]]
    eff = logic.add_efficiency(
        chem.merge(smiles, on="compound_key", how="left"),
        _hac_by_smiles(tuple(sorted(smiles["canonical_smiles"].dropna().unique()))),
    )
    st.plotly_chart(
        charts.efficiency_scatter(
            eff, "ligand_efficiency", "ligand efficiency (kcal/mol per heavy atom)"
        ),
        width="stretch",
    )
    st.divider()
    st.plotly_chart(
        charts.efficiency_scatter(
            eff, "lipophilic_efficiency", "lipophilic efficiency (pChEMBL − logP)"
        ),
        width="stretch",
    )
