# visualization.py

# Standard library imports
import os
from collections import defaultdict
from typing import List, Optional, Union, Tuple, Dict
from importlib import resources

# Third-party imports
import matplotlib.pyplot as plt
import matplotlib.figure
import seaborn as sns
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go


# Internal imports
from .RDDcounts import RDDCounts
from .utils import RDD_counts_to_wide, calculate_proportions


def sort_nodes_by_flow(flows_df, processes_df):
    """
    Sort nodes to minimize crossings in a Sankey diagram.

    Parameters
    ----------
    flows_df : pd.DataFrame
        DataFrame containing 'source', 'target', and 'value' columns.
    processes_df : pd.DataFrame
        DataFrame containing 'id' and 'level' columns for hierarchical levels.

    Returns
    -------
    sorted_nodes : list
        List of node names sorted to minimize crossings.
    node_indices : dict
        Mapping of node names to their sorted indices.
    """
    # Handle missing "id" column by using unique nodes from flows_df
    if "id" not in processes_df.columns:
        if getattr(processes_df.index, "name", None) == "id":
            processes_df = processes_df.reset_index()
        else:
            unique_nodes = pd.concat(
                [flows_df["source"], flows_df["target"]]
            ).unique()
            processes_df = pd.DataFrame(
                {"id": unique_nodes, "level": 0}
            )  # Assign default level

    # Ensure "level" column exists
    if "level" not in processes_df.columns:
        processes_df["level"] = 0  # Default to level 0 if missing

    # Define levels based on the processes DataFrame
    levels_new = processes_df.set_index("id")["level"].to_dict()
    nodes_per_level_new = defaultdict(list)
    for node, level in levels_new.items():
        nodes_per_level_new[level].append(node)

    # Calculate cumulative flow values
    cumulative_outgoing_new = defaultdict(int)
    for _, row in flows_df.iterrows():
        cumulative_outgoing_new[row["source"]] += row["value"]

    # Sort nodes per level
    sorted_nodes_per_level_new = {}

    # First level: Sort by outgoing flow value (descending)
    first_level_new = min(nodes_per_level_new.keys())
    sorted_nodes_per_level_new[first_level_new] = sorted(
        nodes_per_level_new[first_level_new],
        key=lambda x: cumulative_outgoing_new[x],
        reverse=True,
    )

    # Subsequent levels: Align with previous levels
    for level in range(
        first_level_new + 1, max(nodes_per_level_new.keys()) + 1
    ):
        if level in nodes_per_level_new:
            previous_level_nodes = sorted_nodes_per_level_new[level - 1]
            current_level_nodes = nodes_per_level_new[level]

            # Map connections to the current level
            connections_new = defaultdict(list)
            for _, row in flows_df.iterrows():
                if (
                    row["source"] in previous_level_nodes
                    and row["target"] in current_level_nodes
                ):
                    connections_new[row["source"]].append(
                        (row["target"], row["value"])
                    )

            # Sort current level nodes based on previous level
            sorted_current_level_new = []
            for prev_node in previous_level_nodes:
                if prev_node in connections_new:
                    sorted_connections_new = sorted(
                        connections_new[prev_node],
                        key=lambda x: x[1],
                        reverse=True,
                    )
                    sorted_current_level_new.extend(
                        [target for target, _ in sorted_connections_new]
                    )

            # Add unconnected nodes
            remaining_nodes_new = set(current_level_nodes) - set(
                sorted_current_level_new
            )
            sorted_current_level_new.extend(remaining_nodes_new)

            sorted_nodes_per_level_new[level] = sorted_current_level_new

    # Combine all sorted nodes into a single list
    sorted_nodes = []
    for level in sorted(sorted_nodes_per_level_new.keys()):
        sorted_nodes.extend(sorted_nodes_per_level_new[level])

    # Map nodes to indices
    node_indices = {node: idx for idx, node in enumerate(sorted_nodes)}

    return sorted_nodes, node_indices


from abc import ABC, abstractmethod
from typing import Tuple
import pandas as pd


class VisualizationBackend(ABC):  # pragma: no cover
    @abstractmethod
    def plot_reference_type_distribution(
        self,
        data: pd.DataFrame,
        group_by: bool,
        figsize: Tuple[int, int],
        **kwargs,
    ):
        pass

    @abstractmethod
    def box_plot_RDD_proportions(
        self,
        data: pd.DataFrame,
        group_by: bool = False,
        group_colors: Optional[dict] = None,
        figsize: Tuple[int, int] = (10, 6),
        **kwargs,
    ):
        pass

    @abstractmethod
    def plot_RDD_proportion_heatmap(
        self,
        data: pd.DataFrame,
        level: int,
        figsize: Tuple[int, int],
        **kwargs,
    ):
        pass

    @abstractmethod
    def plot_pca_results(
        self,
        pca_df: pd.DataFrame,
        explained_variance: List[float],
        x_pc: str,
        y_pc: str,
        group_by: bool,
        group_colors: Optional[dict],
        figsize: Tuple[int, int],
        **kwargs,
    ):
        pass

    @abstractmethod
    def plot_explained_variance(
        self,
        explained_variance: List[float],
        figsize: Tuple[int, int],
        **kwargs,
    ):
        pass

    @abstractmethod
    def plot_sankey(
        self,
        RDD_counts: "RDDCounts",
        color_mapping_file: str,
        max_hierarchy_level: Optional[int] = None,
        filename_filter: Optional[str] = None,
        dark_mode: bool = False,
    ) -> go.Figure:
        pass


def filter_and_group_RDD_counts(
    RDD_counts_instance: "RDDCounts",
    level: int,
    reference_types: Optional[List[str]] = None,
    group_by: bool = False,
    sample_names: Optional[Union[str, List[str]]] = None,
    top_n: Optional[int] = None,
    top_n_method: str = "per_sample",
    group: Optional[Union[str, List[str]]] = None,
    upper_level: Optional[int] = None,
    lower_level: Optional[int] = None,
    upper_level_reference_types: Optional[List[str]] = None,
) -> pd.DataFrame:
    """
    Filter and group the RDD counts data by ontology level, reference types, and optional parameters.

    Parameters
    ----------
    RDD_counts_instance : RDDCounts
        An instance of the RDDCounts class containing the RDD counts data.
    level : int
        The ontology level to filter by.
    reference_types : list of str, optional
        Specific reference types to include in the filtered data.
    group_by : bool, optional
        Whether to group the data by the 'group' column.
    sample_names : str or list of str, optional
        Filter by specific sample name(s).
    top_n : int, optional
        Select the top N reference types.
    top_n_method : str, optional
        Method to select top N types: 'per_sample', 'total', or 'average'.

    Returns
    -------
    pd.DataFrame
        The filtered and grouped RDD counts data.
    """
    filtered_counts = RDD_counts_instance.filter_counts(
        reference_types=reference_types,
        level=level,
        sample_names=sample_names,
        top_n=top_n,
        top_n_method=top_n_method,
        group=group,
        upper_level=upper_level,
        lower_level=lower_level,
        upper_level_reference_types=upper_level_reference_types,
    )

    if filtered_counts.empty:
        raise ValueError("No data available for the specified filters.")

    if group_by:
        data = (
            filtered_counts.groupby(["reference_type", "group"])["count"]
            .sum()
            .reset_index()
        )
    else:
        data = (
            filtered_counts.groupby("reference_type")["count"]
            .sum()
            .reset_index()
        )

    return data


def prepare_boxplot_data(
    RDD_counts_instance: "RDDCounts",
    level: int,
    reference_types: Optional[List[str]] = None,
    sample_names: Optional[Union[str, List[str]]] = None,
    top_n: Optional[int] = None,
    top_n_method: str = "per_sample",
    group: Optional[Union[str, List[str]]] = None,
    upper_level: Optional[int] = None,
    lower_level: Optional[int] = None,
    upper_level_reference_types: Optional[List[str]] = None,
) -> pd.DataFrame:
    """
    Prepare the RDD proportions data for box plotting.
    """
    # Calculate proportions from complete dataset first
    df_proportions = calculate_proportions(
        RDD_counts_instance.counts, level=level
    )

    # Convert to long format
    df_long = df_proportions.reset_index().melt(
        id_vars=["filename", "group"],
        var_name="reference_type",
        value_name="proportion",
    )

    # Apply hierarchical filtering if specified
    if upper_level is not None and upper_level_reference_types is not None:
        # Get the mapping between upper and lower levels
        ontology_table = RDD_counts_instance.ontology_table
        upper_col = RDD_counts_instance.get_ontology_column_for_level(
            upper_level
        )
        lower_col = RDD_counts_instance.get_ontology_column_for_level(level)

        # Filter to only include reference types under the specified upper level types
        valid_lower_types = ontology_table[
            ontology_table[upper_col].isin(upper_level_reference_types)
        ][lower_col].unique()

        df_long = df_long[df_long["reference_type"].isin(valid_lower_types)]

    # Now filter the long-format proportion data
    if reference_types is not None:
        df_long = df_long[df_long["reference_type"].isin(reference_types)]

    if sample_names is not None:
        if isinstance(sample_names, str):
            sample_names = [sample_names]
        df_long = df_long[df_long["filename"].isin(sample_names)]

    if group is not None:
        if isinstance(group, str):
            group = [group]
        df_long = df_long[df_long["group"].isin(group)]

    # Handle top_n filtering
    if top_n is not None and not df_long.empty:
        if top_n_method == "per_sample":
            # Top N reference types per sample
            top_refs = (
                df_long.sort_values(
                    ["filename", "proportion"], ascending=[True, False]
                )
                .groupby("filename", as_index=False, sort=False)
                .head(top_n)["reference_type"]
                .unique()
            )
        elif top_n_method == "total":
            # Get top N reference types overall
            top_refs = (
                df_long.groupby("reference_type")["proportion"]
                .sum()
                .nlargest(top_n)
                .index.tolist()
            )
        elif top_n_method == "average":
            # Get top N reference types by average proportion
            top_refs = (
                df_long.groupby("reference_type")["proportion"]
                .mean()
                .nlargest(top_n)
                .index.tolist()
            )
        else:
            raise ValueError(f"Unsupported top_n_method: {top_n_method}")

        df_long = df_long[df_long["reference_type"].isin(top_refs)]

    if df_long.empty:
        raise ValueError("No data available after applying filters.")

    return df_long


def prepare_heatmap_data(
    RDD_counts_instance: "RDDCounts",
    level: int,
    reference_types: Optional[List[str]] = None,
) -> pd.DataFrame:
    """
    Prepare the RDD proportions data for the heatmap.

    Parameters
    ----------
    RDD_counts_instance : RDDCounts
        An instance of the RDDCounts class containing the RDD counts data.
    level : int
        The ontology level to filter by.
    reference_types : list of str, optional
        Specific reference types to include in the heatmap. If None, all reference types
        are included. Defaults to None.

    Returns
    -------
    pd.DataFrame
        A DataFrame containing the proportions for the selected reference types and level.
    """
    # Access counts and calculate proportions
    counts = RDD_counts_instance.counts
    df_proportions = calculate_proportions(counts, level=level)

    # Filter by reference types
    if reference_types is not None:
        df_proportions_filtered = df_proportions[
            df_proportions.columns.intersection(reference_types)
        ]
    else:
        df_proportions_filtered = df_proportions

    # Ensure only numeric columns remain
    reference_type_columns = df_proportions_filtered.select_dtypes(
        include=["int", "float"]
    ).columns
    return df_proportions_filtered[reference_type_columns]


class MatplotlibBackend(VisualizationBackend):
    def plot_reference_type_distribution(
        self,
        data: pd.DataFrame,
        group_by: bool,
        figsize: Tuple[int, int],
        **kwargs,
    ):
        fig, ax = plt.subplots(figsize=figsize)

        if group_by:
            sns.barplot(
                x="reference_type",
                y="count",
                hue="group",
                data=data,
                palette="viridis",
                ax=ax,
            )
        else:
            sns.barplot(
                x="reference_type",
                y="count",
                data=data,
                ax=ax,
            )

        ax.set_title("reference type Distribution")
        ax.set_xlabel("reference type")
        ax.set_ylabel("Total Count")
        plt.xticks(rotation=45, ha="right")
        plt.tight_layout()
        return fig

    def box_plot_RDD_proportions(
        self,
        data: pd.DataFrame,
        group_by: bool = False,
        group_colors: Optional[dict] = None,
        figsize: Tuple[int, int] = (10, 6),
        ax=None,
        **kwargs,
    ):
        # Only create figure if ax is not provided
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        else:
            fig = ax.get_figure()

        if group_by:
            # Grouped boxplot using 'hue'
            sns.boxplot(
                x="reference_type",
                y="proportion",
                hue="group",
                data=data,
                ax=ax,
                palette=group_colors,
                showfliers=False,
                orient="v",
                **kwargs,
            )
            sns.stripplot(
                x="reference_type",
                y="proportion",
                hue="group",
                data=data,
                dodge=True,  # Separate points by group
                jitter=True,
                palette=group_colors,
                ax=ax,
                marker="o",
                edgecolor="black",
                alpha=0.7,
                linewidth=0.6,
                legend=False,
            )
        else:
            # Ungrouped boxplot
            sns.boxplot(
                x="reference_type",
                y="proportion",
                data=data,
                ax=ax,
                showfliers=False,
                orient="v",
                **kwargs,
            )
            sns.stripplot(
                x="reference_type",
                y="proportion",
                data=data,
                dodge=False,
                jitter=True,
                ax=ax,
                marker="o",
                edgecolor="black",
                alpha=0.3,
                linewidth=0.6,
            )

        ax.set_title("Proportion Distribution of Selected reference types")
        ax.set_xlabel("reference type")
        ax.set_ylabel("Proportion")
        plt.xticks(rotation=45, ha="right")
        plt.tight_layout()

        return fig

    def plot_RDD_proportion_heatmap(
        self,
        data: pd.DataFrame,
        level: int,
        figsize: Tuple[int, int] = (12, 8),
        **kwargs,
    ):
        """
        Render a heatmap of RDD proportions using Seaborn.

        Parameters
        ----------
        data : pd.DataFrame
            The processed data containing RDD proportions.
        level : int
            The ontology level of the data.
        figsize : tuple of int, optional
            The size of the figure (width, height). Defaults to (12, 8).

        Returns
        -------
        matplotlib.figure.Figure
            The rendered Matplotlib figure.
        """
        plt.figure(figsize=figsize)
        ax = sns.heatmap(
            data, cmap="viridis", annot=False, cbar=True, **kwargs
        )
        ax.set_title(f"Proportion Heatmap of reference types (Level {level})")
        ax.set_xlabel("reference types")
        ax.set_ylabel("Samples")
        plt.xticks(rotation=45, ha="right")
        plt.tight_layout()
        return ax.get_figure()

    def plot_pca_results(
        self,
        pca_df: pd.DataFrame,
        explained_variance: List[float],
        x_pc: str = "PC1",
        y_pc: str = "PC2",
        group_by: bool = True,
        group_colors: Optional[dict] = None,
        figsize: Tuple[int, int] = (10, 6),
        group_column: str = "group",
        **kwargs,
    ):
        """
        Render a PCA results scatter plot using Matplotlib/Seaborn.

        Parameters
        ----------
        pca_df : pd.DataFrame
            DataFrame containing PCA scores and metadata.
        explained_variance : list of float
            List of explained variance ratios for each component.
        x_pc : str, optional
            The principal component to use for the x-axis. Defaults to "PC1".
        y_pc : str, optional
            The principal component to use for the y-axis. Defaults to "PC2".
        group_by : bool, optional
            Whether to color the points by 'group'. Defaults to True.
        group_colors : dict, optional
            A dictionary mapping group names to colors.
        figsize : tuple of int, optional
            The size of the figure (width, height). Defaults to (10, 6).
        group_column : str, optional
            The column name to use for grouping. Defaults to "group".

        Returns
        -------
        matplotlib.figure.Figure
            The rendered Matplotlib figure.
        """
        fig, ax = plt.subplots(figsize=figsize)

        if group_by:
            sns.scatterplot(
                x=x_pc,
                y=y_pc,
                hue=group_column,
                data=pca_df,
                palette=group_colors,
                ax=ax,
            )
        else:
            sns.scatterplot(x=x_pc, y=y_pc, data=pca_df, ax=ax)

        ax.set_title("PCA Plot of RDD Counts")
        ax.set_xlabel(
            f"{x_pc} [{explained_variance[int(x_pc[2]) - 1] * 100:.1f}%]"
        )
        ax.set_ylabel(
            f"{y_pc} [{explained_variance[int(y_pc[2]) - 1] * 100:.1f}%]"
        )
        plt.tight_layout()

        return fig

    def plot_explained_variance(
        self,
        explained_variance: List[float],
        figsize: Tuple[int, int] = (10, 6),
        **kwargs,
    ):
        """
        Render a bar chart for the explained variance using Matplotlib/Seaborn.

        Parameters
        ----------
        explained_variance : list of float
            List of explained variance ratios for each principal component.
        figsize : tuple of int, optional
            The size of the figure (width, height). Defaults to (10, 6).

        Returns
        -------
        matplotlib.figure.Figure
            The rendered Matplotlib figure.
        """
        pc_labels = [f"PC{i+1}" for i in range(len(explained_variance))]
        variance_percentages = [var * 100 for var in explained_variance]

        fig, ax = plt.subplots(figsize=figsize)
        sns.barplot(x=pc_labels, y=variance_percentages, ax=ax, color="blue")

        # Add text annotations
        for i, var in enumerate(variance_percentages):
            ax.text(i, var + 1, f"{var:.1f}%", ha="center")

        ax.set_title("Explained Variance by Principal Component")
        ax.set_xlabel("Principal Component")
        ax.set_ylabel("Explained Variance (%)")
        plt.tight_layout()

        return fig

    def plot_sankey(
        self,
        RDD_counts: "RDDCounts",
        color_mapping_file: str,
        max_hierarchy_level: Optional[int] = None,
        filename_filter: Optional[str] = None,
        dark_mode: bool = False,
    ) -> go.Figure:
        """
        Matplotlib does not support Sankey diagram plotting natively.

        Raises
        ------
        NotImplementedError
            If this method is called.
        """
        raise NotImplementedError(
            "Sankey diagram visualization is not supported in Matplotlib."
        )


class PlotlyBackend(VisualizationBackend):
    def plot_reference_type_distribution(
        self,
        data: pd.DataFrame,
        group_by: bool,
        figsize: Tuple[int, int],
        **kwargs,
    ):
        """
        Render a bar chart for reference type distribution using Plotly.

        Parameters
        ----------
        data : pd.DataFrame
            The filtered and grouped RDD counts data.
        group_by : bool
            Whether to group by the 'group' column.
        figsize : tuple of int
            Ignored for Plotly but maintained for interface compatibility.

        Returns
        -------
        plotly.graph_objects.Figure
            The rendered Plotly figure.
        """
        if group_by:
            fig = px.bar(
                data,
                x="reference_type",
                y="count",
                color="group",
                barmode="group",
                title="reference type Distribution by Group",
                **kwargs,
            )
        else:
            fig = px.bar(
                data,
                x="reference_type",
                y="count",
                title="reference type Distribution",
                **kwargs,
            )

        fig.update_layout(
            xaxis_title="reference type",
            yaxis_title="Total Count",
            xaxis_tickangle=-45,
            template="plotly_white",
        )
        return fig

    def box_plot_RDD_proportions(
        self,
        data: pd.DataFrame,
        group_by: bool = False,
        group_colors: Optional[dict] = None,
        **kwargs,
    ):
        fig = go.Figure()

        if group_by:
            # Iterate over unique groups and add a trace for each
            groups = data["group"].unique()
            for i, group in enumerate(groups):
                group_data = data[data["group"] == group]
                fig.add_trace(
                    go.Box(
                        x=group_data["reference_type"],
                        y=group_data["proportion"],
                        name=group,
                        boxpoints="all",
                        jitter=0.3,
                        pointpos=0,
                        marker=dict(
                            color=(
                                group_colors.get(group)
                                if group_colors
                                else None
                            )
                        ),
                        offsetgroup=i,
                    )
                )
        else:
            # Add one boxplot trace per reference type
            for reference_type in data["reference_type"].unique():
                RDD_data = data[data["reference_type"] == reference_type]
                fig.add_trace(
                    go.Box(
                        x=RDD_data["reference_type"],
                        y=RDD_data["proportion"],
                        name=reference_type,
                        boxpoints="all",
                        jitter=0.3,
                        pointpos=0,
                        marker=dict(
                            color=(
                                group_colors.get(reference_type)
                                if group_colors
                                else None
                            )
                        ),
                    )
                )

        fig.update_layout(
            title="Proportion Distribution of Selected reference types",
            xaxis_title="reference type",
            yaxis_title="Proportion",
            boxmode="group" if group_by else "overlay",
        )
        return fig

    def plot_RDD_proportion_heatmap(
        self,
        data: pd.DataFrame,
        level: int,
        figsize: Tuple[int, int] = (12, 8),
        **kwargs,
    ):
        """
        Render a heatmap of RDD proportions using Plotly.

        Parameters
        ----------
        data : pd.DataFrame
            The processed data containing RDD proportions.
        level : int
            The ontology level of the data.
        figsize : tuple of int, optional
            Ignored for Plotly but maintained for interface compatibility.

        Returns
        -------
        plotly.graph_objects.Figure
            The rendered Plotly figure.
        """
        fig = go.Figure(
            data=go.Heatmap(
                z=data.values,
                x=data.columns,
                y=data.index,
                colorscale="Viridis",
                colorbar=dict(title="Proportion"),
                **kwargs,
            )
        )
        fig.update_layout(
            title=f"Proportion Heatmap of reference types (Level {level})",
            xaxis_title="reference types",
            yaxis_title="Samples",
        )
        return fig

    def plot_pca_results(
        self,
        pca_df: pd.DataFrame,
        explained_variance: List[float],
        x_pc: str = "PC1",
        y_pc: str = "PC2",
        group_by: bool = True,
        group_colors: Optional[dict] = None,
        figsize: Tuple[int, int] = (10, 6),
        group_column: str = "group",
        **kwargs,
    ):
        """
        Render a PCA results scatter plot using Plotly.

        Parameters
        ----------
        pca_df : pd.DataFrame
            DataFrame containing PCA scores and metadata.
        explained_variance : list of float
            List of explained variance ratios for each component.
        x_pc : str, optional
            The principal component to use for the x-axis. Defaults to "PC1".
        y_pc : str, optional
            The principal component to use for the y-axis. Defaults to "PC2".
        group_by : bool, optional
            Whether to color the points by 'group'. Defaults to True.
        group_colors : dict, optional
            A dictionary mapping group names to colors.
        figsize : tuple of int, optional
            Ignored for Plotly but maintained for interface compatibility.
        group_column : str, optional
            The column name to use for grouping. Defaults to "group".

        Returns
        -------
        plotly.graph_objects.Figure
            The rendered Plotly figure.
        """
        fig = go.Figure()

        if group_by:
            groups = pca_df[group_column].unique()
            for group in groups:
                group_data = pca_df[pca_df[group_column] == group]
                fig.add_trace(
                    go.Scatter(
                        x=group_data[x_pc],
                        y=group_data[y_pc],
                        mode="markers",
                        name=group,
                        marker=dict(
                            color=(
                                group_colors.get(group)
                                if group_colors
                                else None
                            )
                        ),
                    )
                )
        else:
            fig.add_trace(
                go.Scatter(
                    x=pca_df[x_pc],
                    y=pca_df[y_pc],
                    mode="markers",
                    name="Samples",
                )
            )

        fig.update_layout(
            title="PCA Plot of RDD Counts",
            xaxis_title=f"{x_pc} [{explained_variance[int(x_pc[2]) - 1] * 100:.1f}%]",
            yaxis_title=f"{y_pc} [{explained_variance[int(y_pc[2]) - 1] * 100:.1f}%]",
            template="plotly_white",
        )

        return fig

    def plot_explained_variance(
        self,
        explained_variance: List[float],
        figsize: Tuple[int, int] = (10, 6),
        **kwargs,
    ):
        """
        Render a bar chart for the explained variance using Plotly.

        Parameters
        ----------
        explained_variance : list of float
            List of explained variance ratios for each principal component.
        figsize : tuple of int, optional
            Ignored for Plotly but maintained for interface compatibility.

        Returns
        -------
        plotly.graph_objects.Figure
            The rendered Plotly figure.
        """
        pc_labels = [f"PC{i+1}" for i in range(len(explained_variance))]
        variance_percentages = [var * 100 for var in explained_variance]

        fig = go.Figure(
            go.Bar(
                x=pc_labels,
                y=variance_percentages,
                text=[f"{var:.1f}%" for var in variance_percentages],
                textposition="auto",
            )
        )
        fig.update_layout(
            title="Explained Variance by Principal Component",
            xaxis_title="Principal Component",
            yaxis_title="Explained Variance (%)",
            template="plotly_white",
        )
        return fig

    def plot_sankey(
        self,
        RDD_counts: "RDDCounts",
        color_mapping_file: str,
        max_hierarchy_level: Optional[int] = None,
        filename_filter: Optional[str] = None,
        dark_mode: bool = False,
    ) -> go.Figure:
        """
        Visualize the RDD flows as a Sankey diagram using Plotly.

        Parameters
        ----------
        RDD_counts : RDDCounts
            An instance of RDDCounts containing the RDD counts data.
        color_mapping_file : str
            CSV file mapping the sample types to their respective colors.
        max_hierarchy_level : int, optional
            Maximum ontology level to calculate flows for. Defaults to all levels.
        filename_filter : str, optional
            Specific sample filename to filter by. If None, uses the full data.
        dark_mode : bool, optional
            If True, apply a dark background theme. Defaults to False.

        Returns
        -------
        plotly.graph_objects.Figure
            The figure object for the Sankey diagram.
        """
        # Generate flows and processes using RDDCounts
        flows_df, processes_df = RDD_counts.generate_RDDflows(
            max_hierarchy_level=max_hierarchy_level,
            filename_filter=filename_filter,
        )

        # Sort nodes to minimize crossings
        sorted_nodes, node_indices = sort_nodes_by_flow(flows_df, processes_df)

        # Map flows to sorted indices
        source_indices = flows_df["source"].map(node_indices)
        target_indices = flows_df["target"].map(node_indices)
        values = flows_df["value"]

        # Load and verify the color mapping
        color_df = pd.read_csv(color_mapping_file, sep=";")
        color_df["color_code"] = color_df["color_code"].fillna("#D3D3D3")
        color_mapping = {
            row["descriptor"]: row["color_code"]
            for _, row in color_df.iterrows()
        }

        # Assign colors to nodes and links
        node_colors = [
            color_mapping.get(node, "#D3D3D3") for node in sorted_nodes
        ]
        link_colors = [
            color_mapping.get(node, "#D3D3D3") for node in flows_df["source"]
        ]

        # Create Sankey diagram
        fig = go.Figure(
            data=[
                go.Sankey(
                    arrangement="fixed",
                    node=dict(
                        pad=15,
                        thickness=20,
                        line=dict(color="black", width=0.5),
                        label=sorted_nodes,
                        color=node_colors,
                    ),
                    link=dict(
                        source=source_indices,
                        target=target_indices,
                        value=values,
                        color=link_colors,
                    ),
                )
            ]
        )

        # Apply dark or light theme
        if dark_mode:
            fig.update_layout(
                title_text="RDD Flows Sankey Diagram",
                font=dict(color="white", size=12),
                paper_bgcolor="black",
                plot_bgcolor="black",
            )
        else:
            fig.update_layout(
                title_text="RDD Flows Sankey Diagram",
                font=dict(color="black", size=12),
                paper_bgcolor="white",
                plot_bgcolor="white",
            )

        return fig


class Visualizer:  # pragma: no cover
    def __init__(self, backend: VisualizationBackend):
        self.backend = backend

    def set_backend(self, backend: VisualizationBackend):
        """
        Set the visualization backend.

        Parameters
        ----------
        backend : VisualizationBackend
            The backend to use for rendering visualizations.
        """
        self.backend = backend

    def plot_reference_type_distribution(
        self,
        RDD_counts_instance: "RDDCounts",
        level: int = 3,
        reference_types: Optional[List[str]] = None,
        group_by: bool = False,
        sample_names: Optional[Union[str, List[str]]] = None,
        top_n: Optional[int] = None,
        top_n_method: str = "per_sample",
        upper_level: Optional[int] = None,
        lower_level: Optional[int] = None,
        upper_level_reference_types: Optional[List[str]] = None,
        group: Optional[Union[str, List[str]]] = None,
        figsize: Tuple[int, int] = (10, 6),
        **kwargs,
    ):
        """
        Plot a bar chart showing the distribution of reference types.

        Parameters
        ----------
        RDD_counts_instance : RDDCounts
            An instance of the RDDCounts class containing the RDD counts data.
        level : int
            The ontology level to filter by.
        reference_types : list of str, optional
            Specific reference types to include in the plot.
        group_by : bool, optional
            Whether to group by the 'group' column in the plot.
        sample_names : str or list of str, optional
            Filter by specific sample name(s).
        top_n : int, optional
            Select the top N reference types.
        top_n_method : str, optional
            Method to select top N reference types.
        figsize : tuple, optional
            Size of the plot figure.

        Returns
        -------
        Figure
            The rendered figure.
        """
        data = filter_and_group_RDD_counts(
            RDD_counts_instance=RDD_counts_instance,
            level=level,
            reference_types=reference_types,
            group_by=group_by,
            sample_names=sample_names,
            top_n=top_n,
            top_n_method=top_n_method,
            group=group,
            upper_level=upper_level,
            lower_level=lower_level,
            upper_level_reference_types=upper_level_reference_types,
        )

        return self.backend.plot_reference_type_distribution(
            data, group_by=group_by, figsize=figsize, **kwargs
        )

    def box_plot_RDD_proportions(
        self,
        RDD_counts_instance: "RDDCounts",
        level: int = 3,
        reference_types: Optional[List[str]] = None,
        group_by: bool = False,
        group_colors: Optional[dict] = None,
        figsize: Tuple[int, int] = (10, 6),
        sample_names: Optional[Union[str, List[str]]] = None,
        top_n: Optional[int] = None,
        top_n_method: str = "per_sample",
        group: Optional[Union[str, List[str]]] = None,
        upper_level: Optional[int] = None,
        lower_level: Optional[int] = None,
        upper_level_reference_types: Optional[List[str]] = None,
        **kwargs,
    ):
        """
        Plot box plots showing the distribution of RDD proportions.

        Parameters
        ----------
        RDD_counts_instance : RDDCounts
            An instance of the RDDCounts class containing the RDD counts data.
        level : int, optional
            The ontology level to filter by. Defaults to 3.
        reference_types : list of str, optional
            Specific reference types to include. Defaults to None.
        group_by : bool, optional
            If True, groups the data by the 'group' column. Defaults to False.
        group_colors : dict, optional
            A dictionary mapping group names to colors. Defaults to None.
        figsize : tuple of int, optional
            The size of the figure (width, height). Defaults to (10, 6).
        sample_names : str or list of str, optional
            Filter by specific sample name(s).
        top_n : int, optional
            Select top N reference types.
        top_n_method : str, optional
            Selection method: 'per_sample', 'total', or 'average'.

        Returns
        -------
        matplotlib.figure.Figure or plotly.graph_objects.Figure
            The rendered figure object.
        """
        data = prepare_boxplot_data(
            RDD_counts_instance=RDD_counts_instance,
            level=level,
            reference_types=reference_types,
            sample_names=sample_names,
            top_n=top_n,
            top_n_method=top_n_method,
            group=group,
            upper_level=upper_level,
            lower_level=lower_level,
            upper_level_reference_types=upper_level_reference_types,
        )

        return self.backend.box_plot_RDD_proportions(
            data,
            group_by=group_by,
            group_colors=group_colors,
            figsize=figsize,
            **kwargs,
        )

    def plot_RDD_proportion_heatmap(
        self,
        RDD_counts_instance: "RDDCounts",
        level: int = 3,
        reference_types: Optional[List[str]] = None,
        figsize: Tuple[int, int] = (12, 8),
        **kwargs,
    ):
        """
        Plot a heatmap of RDD proportions for selected reference types.

        Parameters
        ----------
        RDD_counts_instance : RDDCounts
            An instance of the RDDCounts class containing the RDD counts data.
        level : int, optional
            The ontology level to filter by. Defaults to 3.
        reference_types : list of str, optional
            Specific reference types to include in the heatmap. Defaults to None.
        figsize : tuple of int, optional
            The size of the figure (width, height). Defaults to (12, 8).

        Returns
        -------
        matplotlib.figure.Figure or plotly.graph_objects.Figure
            The rendered heatmap.
        """
        # Prepare data
        data = prepare_heatmap_data(
            RDD_counts_instance, level, reference_types
        )

        # Render using the backend
        return self.backend.plot_RDD_proportion_heatmap(
            data, level=level, figsize=figsize, **kwargs
        )

    def plot_pca_results(
        self,
        pca_df: pd.DataFrame,
        explained_variance: List[float],
        x_pc: str = "PC1",
        y_pc: str = "PC2",
        group_by: bool = True,
        group_colors: Optional[dict] = None,
        figsize: Tuple[int, int] = (10, 6),
        group_column: str = "group",
        **kwargs,
    ):
        """
        Plot the PCA results using the selected backend.

        Parameters
        ----------
        pca_df : pd.DataFrame
            DataFrame containing PCA scores and metadata.
        explained_variance : list of float
            List of explained variance ratios for each component.
        x_pc : str, optional
            The principal component to use for the x-axis. Defaults to "PC1".
        y_pc : str, optional
            The principal component to use for the y-axis. Defaults to "PC2".
        group_by : bool, optional
            Whether to color the points by 'group'. Defaults to True.
        group_colors : dict, optional
            A dictionary mapping group names to colors.
        figsize : tuple of int, optional
            The size of the figure (width, height). Defaults to (10, 6).
        group_column : str, optional
            The column name to use for grouping. Defaults to "group".

        Returns
        -------
        matplotlib.figure.Figure or plotly.graph_objects.Figure
            The rendered PCA plot.
        """
        return self.backend.plot_pca_results(
            pca_df=pca_df,
            explained_variance=explained_variance,
            x_pc=x_pc,
            y_pc=y_pc,
            group_by=group_by,
            group_colors=group_colors,
            figsize=figsize,
            group_column=group_column,
            **kwargs,
        )

    def plot_explained_variance(
        self,
        explained_variance: List[float],
        figsize: Tuple[int, int] = (10, 6),
        **kwargs,
    ):
        """
        Plot a bar chart of explained variance for principal components.

        Parameters
        ----------
        explained_variance : list of float
            List of explained variance ratios for each principal component.
        figsize : tuple of int, optional
            The size of the figure (width, height). Defaults to (10, 6).

        Returns
        -------
        matplotlib.figure.Figure or plotly.graph_objects.Figure
            The rendered bar chart.
        """
        return self.backend.plot_explained_variance(
            explained_variance=explained_variance, figsize=figsize, **kwargs
        )

    def plot_sankey(
        self,
        RDD_counts: "RDDCounts",
        color_mapping_file: str,
        max_hierarchy_level: Optional[int] = None,
        filename_filter: Optional[str] = None,
        dark_mode: bool = False,
    ) -> go.Figure:
        """
        Plot the Sankey diagram using the selected backend.

        Parameters
        ----------
        RDD_counts : RDDCounts
            An instance of RDDCounts containing the RDD counts data.
        color_mapping_file : str
            CSV file mapping the sample types to their respective colors.
        max_hierarchy_level : int, optional
            Maximum ontology level to calculate flows for. Defaults to all levels.
        filename_filter : str, optional
            Specific sample filename to filter by. If None, uses the full data.
        dark_mode : bool, optional
            If True, apply a dark background theme. Defaults to False.

        Returns
        -------
        plotly.graph_objects.Figure
            The figure object for the Sankey diagram.
        """
        return self.backend.plot_sankey(
            RDD_counts=RDD_counts,
            color_mapping_file=color_mapping_file,
            max_hierarchy_level=max_hierarchy_level,
            filename_filter=filename_filter,
            dark_mode=dark_mode,
        )
