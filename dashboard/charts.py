"""Plotly figure builders for the SARVault dashboard (pure: DataFrame -> Figure)."""

import plotly.express as px
import plotly.graph_objects as go


def sar_ranking_bar(df, top_n: int = 20) -> go.Figure:
    """Horizontal bar of the top-N compounds by median pchembl."""
    top = df.sort_values("median_pchembl", ascending=False).head(top_n)
    fig = px.bar(
        top,
        x="median_pchembl",
        y="molecule_chembl_id",
        orientation="h",
        hover_data=["max_pchembl", "n_measurements", "n_assays"],
        labels={"median_pchembl": "median pChEMBL", "molecule_chembl_id": "compound"},
    )
    fig.update_layout(yaxis={"categoryorder": "total ascending"}, height=520)
    return fig


def selectivity_scatter(df) -> go.Figure:
    """Selectivity index vs. best-target potency for multi-target compounds."""
    fig = px.scatter(
        df,
        x="best_pchembl",
        y="selectivity_index",
        hover_data=["molecule_chembl_id", "best_target", "n_targets"],
        labels={
            "best_pchembl": "best-target median pChEMBL",
            "selectivity_index": "selectivity (log10 fold)",
        },
    )
    fig.update_layout(height=520)
    return fig


def chemical_space_scatter(df, x: str = "mw_freebase", y: str = "alogp") -> go.Figure:
    """Chemical-space scatter, coloured by approval and sized by potency."""
    plot_df = df.dropna(subset=[x, y, "best_pchembl"])
    fig = px.scatter(
        plot_df,
        x=x,
        y=y,
        color="is_approved_drug",
        size="best_pchembl",
        hover_data=["molecule_chembl_id", "pref_name"],
    )
    fig.update_layout(height=520)
    return fig


def efficiency_scatter(df, metric: str, label: str) -> go.Figure:
    """Ligand-efficiency metric vs. best-target potency, coloured by approval.

    ``metric`` is a column such as 'ligand_efficiency' or 'lipophilic_efficiency';
    points to the upper right are both potent and efficient.
    """
    plot_df = df.dropna(subset=["best_pchembl", metric])
    fig = px.scatter(
        plot_df,
        x="best_pchembl",
        y=metric,
        color="is_approved_drug",
        hover_data=["molecule_chembl_id", "heavy_atoms", "mw_freebase"],
        labels={"best_pchembl": "best-target median pChEMBL", metric: label},
    )
    fig.update_layout(height=440)
    return fig


def property_histogram(df, column: str = "mw_freebase") -> go.Figure:
    """Overlaid histogram of a physicochemical property by approval status."""
    fig = px.histogram(df, x=column, color="is_approved_drug", barmode="overlay", nbins=40)
    fig.update_layout(height=360)
    return fig


def sar_heatmap(df, top_n: int = 30) -> go.Figure:
    """Compound × target median-pChEMBL heatmap for the top-N compounds.

    Compounds are ranked by their best (max) median potency across the in-scope
    targets. The pivot leaves gaps where a compound was not measured against a
    target; those read as blank cells. One bright cell = selective, several =
    promiscuous.
    """
    if df.empty:
        return go.Figure()
    ranked = df.groupby("molecule_chembl_id")["median_pchembl"].max().sort_values(ascending=False)
    top_ids = ranked.head(top_n).index.tolist()
    wide = (
        df[df["molecule_chembl_id"].isin(top_ids)]
        .pivot_table(
            index="molecule_chembl_id",
            columns="target_pref_name",
            values="median_pchembl",
            aggfunc="max",
        )
        .reindex(top_ids)
    )
    fig = px.imshow(
        wide,
        color_continuous_scale="Viridis",
        aspect="auto",
        labels={"x": "target", "y": "compound", "color": "median pChEMBL"},
    )
    fig.update_layout(
        height=max(360, 22 * len(top_ids) + 120),
        margin={"l": 10, "r": 10, "t": 10, "b": 40},
    )
    return fig


def target_potency_violin(df) -> go.Figure:
    """Distribution (violin + inner box) of per-compound median pChEMBL per target."""
    fig = px.violin(
        df,
        x="target_pref_name",
        y="median_pchembl",
        box=True,
        points=False,
        labels={"target_pref_name": "target", "median_pchembl": "median pChEMBL"},
    )
    fig.update_layout(height=420, margin={"l": 10, "r": 10, "t": 10, "b": 40})
    return fig


def compound_potency_bar(df) -> go.Figure:
    """Per-target median potency for a single compound; height scales with targets."""
    fig = px.bar(
        df,
        x="median_pchembl",
        y="target",
        orientation="h",
        hover_data=["max_pchembl", "n_measurements"],
        labels={"median_pchembl": "median pChEMBL", "target": ""},
    )
    fig.update_layout(
        yaxis={"categoryorder": "total ascending"},
        height=max(140, 46 * len(df) + 70),
        margin={"l": 10, "r": 10, "t": 10, "b": 40},
    )
    return fig


def standard_type_bar(df) -> go.Figure:
    """Bar chart of activities per ChEMBL standard (endpoint) type: IC50, Ki, ..."""
    fig = px.bar(
        df,
        x="standard_type",
        y="n_activities",
        hover_data=["n_compounds"],
        labels={"standard_type": "standard type", "n_activities": "activities"},
    )
    fig.update_layout(
        height=360, xaxis={"type": "category", "categoryorder": "total descending"}
    )
    return fig


def confidence_bar(df) -> go.Figure:
    """Bar chart of activities per ChEMBL assay confidence score."""
    fig = px.bar(
        df,
        x="confidence_score",
        y="n_activities",
        hover_data=["n_assays"],
        labels={"confidence_score": "confidence score", "n_activities": "activities"},
    )
    fig.update_layout(height=360, xaxis={"type": "category"})
    return fig
