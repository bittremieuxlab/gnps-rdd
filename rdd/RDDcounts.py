# RDDcounts.py

# Standard library imports
import os
from typing import List, Optional, Tuple, Union

# Third-party imports
import numpy as np
import pandas as pd

# Internal imports
from .utils import (
    _load_RDD_metadata,
    _load_sample_types,
    _validate_groups,
    get_sample_metadata,
    normalize_network,
    split_reference_sample,
    get_gnps_task_data,
)


class RDDCounts:
    """
    A class for generating Reference Data-Driven (RDD) counts across ontology levels
    using GNPS molecular networking data and reference metadata.

    Parameters
    ----------
    gnps_network_path : str
        Path to the GNPS network file (.tsv format).
    sample_types : str
        Indicates whether to use 'all', 'simple', or 'complex' references.
    sample_groups : list of str, optional
        Sample group identifiers to include (e.g., ["G1", "G2"]).
    reference_groups : list of str, optional
        Reference group identifiers to include (e.g., ["G3"]).
    sample_group_col : str
        The column in the sample metadata representing sample group identifiers.
        Default is "group".
    levels : int, optional
        Number of ontology levels to analyze. Default is 6.
        If `ontology_columns` is provided, this still sets how many levels are analyzed.
    external_reference_metadata : str, optional
        Path to a user-supplied reference metadata file.
        If None, internal foodomics metadata is used.
    external_sample_metadata : str, optional
        Path to a user-supplied sample metadata file. If None, metadata is extracted from GNPS.
    ontology_columns : list of str, optional
        Ordered list of ontology column names to use instead of default `sample_type_groupX` columns.
        These will be renamed internally to `<column_name><level_number>` to support level-based logic.
    """

    def __init__(
        self,
        sample_types: str,
        gnps_network_path: Optional[str] = None,
        sample_groups: Optional[List[str]] = None,
        reference_groups: Optional[List[str]] = None,
        sample_group_col: str = "group",
        levels: int = 6,
        external_reference_metadata: Optional[str] = None,
        external_sample_metadata: Optional[str] = None,
        ontology_columns: Optional[List[str]] = None,
        blank_identifier: str = "water",
        task_id: Optional[str] = None,
        gnps_2: bool = True,
    ) -> None:

        if (task_id is None) == (gnps_network_path is None):
            raise ValueError(
                "Provide exactly one of task_id or gnps_network_path."
            )
        if task_id is None:
            self.raw_gnps_network = pd.read_csv(gnps_network_path, sep="\t")
        else:
            gnps_data = get_gnps_task_data(task_id, gnps_2)
            self.raw_gnps_network = gnps_data
        self.sample_types = sample_types
        self.sample_groups = sample_groups
        self.reference_groups = reference_groups
        self.levels = levels
        self.sample_group_col = sample_group_col
        self.blank_identifier = blank_identifier

        self.reference_metadata = _load_RDD_metadata(
            external_reference_metadata
        )
        self.sample_metadata = get_sample_metadata(
            self.raw_gnps_network,
            sample_groups,
            external_sample_metadata,
            filename_col="filename",
        )
        self.sample_types_df, self.ontology_columns_renamed = (
            _load_sample_types(
                self.reference_metadata,
                sample_types,
                ontology_columns=ontology_columns,
            )
        )

        if self.ontology_columns_renamed:
            if self.levels > len(self.ontology_columns_renamed):
                raise ValueError(
                    f"levels ({self.levels}) exceeds provided ontology columns "
                    f"({len(self.ontology_columns_renamed)})."
                )

        self.ontology_table = (
            self.sample_types_df.copy()
            .drop(columns=["sample_name"], errors="ignore")
            .drop_duplicates()
        )
        self.normalized_network = normalize_network(
            self.raw_gnps_network, self.sample_groups, self.reference_groups
        )
        self.file_level_counts = self.file_counts(
            self.normalized_network,
            self.reference_metadata,
            self.sample_metadata,
        )
        self.counts = self.create_RDD_counts_all_levels()

    def _get_ontology_column_for_level(self, level: int) -> str:
        """
        Retrieves the appropriate ontology column name corresponding to a given level.

        This method returns the renamed ontology column (if custom ontology columns were
        provided and renamed), or defaults to the standard column format
        'sample_type_group{level}'.

        Parameters
        ----------
        level : int
            The ontology level for which to retrieve the column name.

        Returns
        -------
        str
            The name of the ontology column corresponding to the specified level.
        """
        if self.ontology_columns_renamed:
            if level < 1 or level > len(self.ontology_columns_renamed):
                raise ValueError(
                    f"Invalid level {level}; expected 1..{len(self.ontology_columns_renamed)}."
                )
            return self.ontology_columns_renamed[level - 1]
        return f"sample_type_group{level}"

    def file_counts(
        self,
        normalized_gnps_network,
        reference_metadata,
        sample_metadata,
        reference_name_col="sample_name",
    ):
        """
        Compute file-level RDD counts by matching sample and reference files to shared molecular clusters.

        This function identifies shared clusters between sample and reference files,
        counts the number of shared clusters per (filename, reference_type) pair,
        and returns the data in long format. It assigns level 0 and attaches group labels
        from the sample metadata.

        Parameters
        ----------
        normalized_gnps_network : pd.DataFrame
            Pre-processed GNPS network in long format with 'filename' and 'cluster_index'.
        reference_metadata : pd.DataFrame
            Reference metadata containing at least ['filename', reference_name_col].
        sample_metadata : pd.DataFrame
            Sample metadata containing at least ['filename', sample_group_col].
            The column in `sample_metadata` representing sample group identifiers. Default is "group".
        reference_name_col : str, optional
            The column in `reference_metadata` that serves as the reference name (e.g., 'sample_name'). Default is "sample_name".

        Returns
        -------
        pd.DataFrame
            A long-format DataFrame containing file-level RDD counts with columns:
            - filename : str
            - reference_type : str
            - count : int
            - level : int (always 0 for file-level)
            - group : str
        """

        sample_clusters, reference_clusters = split_reference_sample(
            normalized_gnps_network,
            reference_metadata,
            sample_metadata,
            sample_group_col=self.sample_group_col,
            reference_name_col=reference_name_col,
        )
        sample_clusters.drop_duplicates(
            subset=["filename", "cluster_index"], inplace=True
        )
        reference_clusters.drop_duplicates(
            subset=["filename", reference_name_col, "cluster_index"],
            inplace=True,
        )
        shared_clusters = reference_clusters.merge(
            sample_clusters,
            on="cluster_index",
            suffixes=("_reference", "_sample"),
        )
        # Remove duplicates to avoid double counting
        shared_clusters.drop_duplicates(
            subset=["filename_reference", "filename_sample", "cluster_index"],
            inplace=True,
        )

        self.shared_clusters = shared_clusters

        cluster_count = (
            shared_clusters.groupby(["filename_sample", reference_name_col])
            .size()
            .unstack(fill_value=0)
        )
        cluster_count = cluster_count.drop_duplicates()
        cluster_count_long = cluster_count.stack().reset_index(name="count")
        cluster_count_long.rename(
            columns={
                reference_name_col: "reference_type",
                "filename_sample": "filename",
            },
            inplace=True,
        )
        cluster_count_long["level"] = 0
        filename_to_group = self.sample_metadata.set_index("filename")[
            self.sample_group_col
        ].to_dict()
        cluster_count_long["group"] = cluster_count_long["filename"].map(
            filename_to_group
        )
        return cluster_count_long

    def create_RDD_counts_all_levels(self) -> pd.DataFrame:
        """
        Generate RDD counts across all ontology levels, comparing to 'water' and returning filtered results.

        For each ontology level, this function:
        - Maps `reference_type` entries to their respective ontology groups.
        - Filters out counts less frequent than 'water'.
        - Returns a long-format DataFrame with counts for each (filename, reference_type, level).

        Returns
        -------
        pd.DataFrame
            A DataFrame containing RDD counts with the following columns:
            - filename : str
            - reference_type : str
            - count : int
            - level : int
            - group : str
        """
        RDD_counts_file_level = self.file_level_counts
        RDD_counts_all_levels = [
            RDD_counts_file_level
        ]  # Initialize a list for storing data at all levels
        if "reference_type" not in RDD_counts_file_level.columns:
            raise ValueError(
                "Expected 'reference_type' column in file-level counts."
            )
        RDD_counts_file_level_sample_types = RDD_counts_file_level.merge(
            self.sample_types_df,
            left_on="reference_type",
            right_on="sample_name",
        ).drop_duplicates()

        sample_metadata_map = self.sample_metadata.set_index("filename")[
            self.sample_group_col
        ].to_dict()  # Create the mapping once

        for level in range(1, self.levels + 1):

            ontology_col = self._get_ontology_column_for_level(level)

            RDD_counts_level = (
                RDD_counts_file_level_sample_types.groupby(
                    ["filename", ontology_col]
                )["count"]
                .sum()
                .reset_index()
            )

            wide_format_counts = RDD_counts_level.pivot_table(
                index="filename",
                columns=ontology_col,
                values="count",
                fill_value=0,
            ).reset_index()

            if self.blank_identifier in wide_format_counts.columns:
                water_counts = wide_format_counts[self.blank_identifier]
                columns_to_modify = wide_format_counts.columns.difference(
                    ["filename", self.blank_identifier]
                )
                wide_format_counts.loc[
                    :, columns_to_modify
                ] = wide_format_counts.loc[:, columns_to_modify].where(
                    wide_format_counts.loc[:, columns_to_modify].gt(
                        water_counts, axis=0
                    ),
                    0,
                )
                wide_format_counts = wide_format_counts.drop(
                    columns=[self.blank_identifier]
                )

            wide_format_counts = wide_format_counts.loc[
                :, (wide_format_counts != 0).any(axis=0)
            ]
            if wide_format_counts.empty:
                continue  # Skip this level
            RDD_counts_level = wide_format_counts.melt(
                id_vars="filename",
                var_name="reference_type",
                value_name="count",
            )
            RDD_counts_level["level"] = level

            RDD_counts_all_levels.append(RDD_counts_level)

        RDD_counts_all_levels = pd.concat(
            RDD_counts_all_levels, ignore_index=True
        )

        # Map group information from the sample_metadata to the final DataFrame
        RDD_counts_all_levels["group"] = RDD_counts_all_levels["filename"].map(
            sample_metadata_map
        )

        # Cast 'count' as an integer
        RDD_counts_all_levels["count"] = RDD_counts_all_levels["count"].astype(
            int
        )

        return RDD_counts_all_levels

    def filter_counts(
        self,
        reference_types: Optional[List[str]] = None,
        level: int = 3,
        sample_names: Optional[Union[str, List[str]]] = None,
        group: Optional[Union[str, List[str]]] = None,
        top_n: Optional[int] = None,
        top_n_method: str = "per_sample",
        upper_level: Optional[int] = None,
        lower_level: Optional[int] = None,
        upper_level_reference_types: Optional[List[str]] = None,
    ) -> pd.DataFrame:
        """
        Filters the RDD counts by reference types, ontology level, sample names, groups,
        and optionally by top N types.

        Parameters
        ----------
        reference_types : list of str, optional
            List of reference types to filter by.
        level : int, optional
            Ontology level to filter by. Defaults to 3.
        sample_names : str or list of str, optional
            Filter by specific sample name(s).
        group : str or list of str, optional
            Filter by specific sample group(s).
        top_n : int, optional
            Select the top N reference types based on the method defined in `top_n_method`.
        top_n_method : str, optional
            Method to select top N reference types. Options:
            - 'per_sample': Top N per sample (union across samples)
            - 'total': Top N by total counts across all samples
            - 'average': Top N by mean count per sample

        Returns
        -------
        pd.DataFrame
            Filtered RDD counts dataframe.
        """
        if self.counts is None:
            raise ValueError("RDD counts are not initialized.")

        if upper_level is not None and lower_level is not None:
            if upper_level >= lower_level:
                raise ValueError("upper_level must be lower than lower_level.")

            if upper_level_reference_types is None:
                raise ValueError(
                    "Must provide upper_level_reference_types when using upper_level filtering."
                )

            # Resolve level column names via helper to support both default and custom columns
            upper_ontology_col = self._get_ontology_column_for_level(
                upper_level
            )
            lower_ontology_col = self._get_ontology_column_for_level(
                lower_level
            )

            # Matching lower-level reference types
            reference_types = (
                self.ontology_table[
                    self.ontology_table[upper_ontology_col].isin(
                        upper_level_reference_types
                    )
                ][lower_ontology_col]
                .dropna()
                .unique()
                .tolist()
            )

            # Override level to lower_level for filtering
            level = lower_level

        # Filter by ontology level
        filtered_df = self.counts[self.counts["level"] == level]

        # Filter by sample name(s)
        if sample_names:
            if isinstance(sample_names, str):
                sample_names = [sample_names]
            filtered_df = filtered_df[
                filtered_df["filename"].isin(sample_names)
            ]

        # Filter by group(s)
        if group is not None:
            if isinstance(group, str):
                group = [group]
            filtered_df = filtered_df[filtered_df["group"].isin(group)]

        # Select top N reference types if requested
        if top_n is not None:
            if top_n_method == "per_sample":
                top_df = (
                    filtered_df.sort_values(
                        ["filename", "count"], ascending=[True, False]
                    )
                    .groupby("filename")
                    .head(top_n)
                    .reset_index(drop=True)
                )
                top_reference_types = (
                    top_df["reference_type"].dropna().unique().tolist()
                )

            elif top_n_method == "total":
                top_reference_types = (
                    filtered_df.groupby("reference_type")["count"]
                    .sum()
                    .nlargest(top_n)
                    .index.tolist()
                )

            elif top_n_method == "average":
                top_reference_types = (
                    filtered_df.groupby("reference_type")["count"]
                    .mean()
                    .nlargest(top_n)
                    .index.tolist()
                )

            else:
                raise ValueError(
                    "Invalid top_n_method. Choose from 'per_sample', 'total', or 'average'."
                )

            filtered_df = filtered_df[
                filtered_df["reference_type"].isin(top_reference_types)
            ]

        # Filter again by explicitly provided reference_types
        if reference_types is not None:
            filtered_df = filtered_df[
                filtered_df["reference_type"].isin(reference_types)
            ]

        return filtered_df

    def update_groups(
        self,
        metadata_source: Union[str, dict],
        merge_column: Optional[str] = None,
    ) -> None:
        """
        Updates the 'group' column in the RDD counts and sample_metadata
        DataFrames based on user-provided metadata.

        Parameters
        ----------
        metadata_source : str or dict
            Either:
            - Path to the metadata file (CSV or TSV) containing updated group information, or
            - Dictionary mapping from current group values to new group values
            (e.g., {"G1": "Vegan", "G2": "Omnivore"})
        merge_column : str, optional
            The column in the metadata file to use for updating the group information.
            Required when metadata_source is a file path, ignored when it's a dictionary.

        Raises
        ------
        ValueError
            If the metadata file is not a valid CSV or TSV file, if the necessary
            columns are missing, or if merge_column is not provided for file input.
        """

        if isinstance(metadata_source, dict):
            # Dictionary mapping: map current group values to new ones
            group_mapping = metadata_source

            # Update the 'group' column in counts using the mapping
            self.counts["group"] = (
                self.counts["group"]
                .map(group_mapping)
                .fillna(self.counts["group"])
            )

            # Update the sample_metadata using the sample_group_col
            self.sample_metadata[self.sample_group_col] = (
                self.sample_metadata[self.sample_group_col]
                .map(group_mapping)
                .fillna(self.sample_metadata[self.sample_group_col])
            )

        else:
            # File-based mapping: existing functionality
            if merge_column is None:
                raise ValueError(
                    "merge_column must be provided when metadata_source is a file path."
                )

            metadata_file = metadata_source

            # Detect file extension to determine the separator
            file_extension = os.path.splitext(metadata_file)[1].lower()
            if file_extension == ".csv":
                metadata = pd.read_csv(metadata_file)
            elif file_extension in [".tsv", ".txt"]:
                metadata = pd.read_csv(metadata_file, sep="\t")
            else:
                raise ValueError("Metadata file must be either a CSV or TSV.")

            # Ensure that the metadata contains the necessary columns
            if not {"filename", merge_column}.issubset(metadata.columns):
                raise ValueError(
                    f"Metadata file must contain 'filename' and '{merge_column}' columns."
                )

            # Create a mapping from filename to new group
            filename_to_group = metadata.set_index("filename")[
                merge_column
            ].to_dict()

            # Update the 'group' column in counts
            self.counts["group"] = (
                self.counts["filename"]
                .map(filename_to_group)
                .fillna(self.counts["group"])
            )

            # Update the sample_metadata
            self.sample_metadata[self.sample_group_col] = (
                self.sample_metadata["filename"]
                .map(filename_to_group)
                .fillna(self.sample_metadata[self.sample_group_col])
            )

    def generate_RDDflows(
        self,
        max_hierarchy_level: Optional[int] = None,
        filename_filter: Optional[str] = None,
    ) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Generates RDD flows and processes across ontology levels, with options for
        filtering by a specific sample and specifying a maximum hierarchy level.

        Parameters
        ----------
        max_hierarchy_level : int, optional
            The maximum level to calculate flows up to. Defaults to the instance's levels.
        filename_filter : str, optional
            A specific sample filename to filter by. If None, the full counts dataframe is used.

        Returns
        -------
        Tuple[pd.DataFrame, pd.DataFrame]
            A tuple containing:
            - flows: DataFrame with source, target, and value columns.
            - processes: DataFrame with unique nodes across the levels.
        """
        # Use provided max_hierarchy_level or default to the instance's levels
        max_level = (
            max_hierarchy_level
            if max_hierarchy_level is not None
            else self.levels
        )

        # Filter counts by filename if a filter is specified
        if filename_filter:
            counts = self.counts[self.counts["filename"] == filename_filter]
        else:
            counts = self.counts

        flows = []

        for i in range(1, max_level):
            source_level = i
            target_level = i + 1

            # Dynamically determine column names based on renaming
            if self.ontology_columns_renamed:
                source_col = self.ontology_columns_renamed[source_level - 1]
                target_col = self.ontology_columns_renamed[target_level - 1]
            else:
                source_col = f"sample_type_group{source_level}"
                target_col = f"sample_type_group{target_level}"

            # Group RDD counts at target level
            target_counts = (
                counts[counts["level"] == target_level]
                .groupby("reference_type")["count"]
                .sum()
                .reset_index()
            )

            # Merge with sample types to get source-level info
            merged_df = pd.merge(
                target_counts,
                self.sample_types_df,
                left_on="reference_type",
                right_on=target_col,
            )[[source_col, "reference_type", "count"]].drop_duplicates()

            # Format source-target pairs with levels
            flow = merged_df.rename(
                columns={
                    source_col: "source",
                    "reference_type": "target",
                }
            )
            flow["source"] = flow["source"] + f"_{source_level}"
            flow["target"] = flow["target"] + f"_{target_level}"
            flow.rename(columns={"count": "value"}, inplace=True)

            flows.append(flow)

        # Concatenate flows into a single DataFrame
        flows_df = pd.concat(flows, ignore_index=True)

        # Build processes from unique nodes in flows
        all_nodes = (
            pd.concat([flows_df["source"], flows_df["target"]])
            .dropna()
            .unique()
        )
        processes_df = pd.DataFrame(
            {
                "id": all_nodes,
                "level": [int(node.split("_")[-1]) for node in all_nodes],
            }
        ).set_index("id")

        return flows_df, processes_df
