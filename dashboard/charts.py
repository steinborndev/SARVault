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


def property_histogram(df, column: str = "mw_freebase") -> go.Figure:
    """Overlaid histogram of a physicochemical property by approval status."""
    fig = px.histogram(df, x=column, color="is_approved_drug", barmode="overlay", nbins=40)
    fig.update_layout(height=360)
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
