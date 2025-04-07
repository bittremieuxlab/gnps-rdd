# RDDcounts.py

# Standard library imports
import os
from typing import List, Optional, Tuple

# Third-party imports
import numpy as np
import pandas as pd

# Internal imports
from utils import (
    _load_RDD_metadata,
    _load_sample_types,
    _validate_groups,
    get_sample_metadata,
    normalize_network,
    split_reference_sample,
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
        gnps_network_path: str,
        sample_types: str,
        sample_groups: List[str] = None,
        reference_groups: List[str] = None,
        sample_group_col: str = "group",
        levels: int = 6,
        external_reference_metadata: Optional[str] = None,
        external_sample_metadata: Optional[str] = None,
        ontology_columns: Optional[List[str]] = None,
    ) -> None:

        self.raw_gnps_network = pd.read_csv(gnps_network_path, sep="\t")
        self.sample_types = sample_types
        self.sample_groups = sample_groups
        self.reference_groups = reference_groups
        self.levels = levels
        self.sample_group_col = sample_group_col

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
        cluster_count_long["group"] = cluster_count_long.merge(
            self.sample_metadata, on="filename", how="inner"
        )[self.sample_group_col]
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

            if "water" in wide_format_counts.columns:
                water_counts = wide_format_counts["water"]
                columns_to_modify = wide_format_counts.columns.difference(
                    ["filename", "water"]
                )
                wide_format_counts.loc[
                    :, columns_to_modify
                ] = wide_format_counts.loc[:, columns_to_modify].where(
                    wide_format_counts.loc[:, columns_to_modify].gt(
                        water_counts, axis=0
                    ),
                    0,
                )
                wide_format_counts = wide_format_counts.drop(columns=["water"])

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
        self, reference_types: Optional[List[str]] = None, level: int = 3
    ) -> pd.DataFrame:
        """
        Filters the RDD counts based on reference types and ontology level.

        Parameters
        ----------
        reference_types : list of str, optional
            List of reference types to filter by. If None, all reference types at the specified
            level are included. Defaults to None.
        level : int, optional
            The ontology level to filter the RDD counts by. Defaults to 3.

        Returns
        -------
        pd.DataFrame
            A DataFrame containing the filtered RDD counts with columns: filename,
            reference type, level, count, and group.

        Raises
        ------
        ValueError
            If the  RDD counts have not been yet created
        """
        if self.counts is None:
            raise ValueError(
                "RDD counts have not been created yet. Call create() first."
            )
        if reference_types is None:
            filtered_df = self.counts[self.counts["level"] == level]
            return filtered_df
        filtered_df = self.counts[
            (self.counts["reference_type"].isin(reference_types))
            & (self.counts["level"] == level)
        ]
        return filtered_df

    def update_groups(self, metadata_file: str, merge_column: str) -> None:
        """
        Updates the 'group' column in the RDD counts and sample_metadata
        DataFrames based on user-provided metadata.

        Parameters
        ----------
        metadata_file : str
            Path to the metadata file (CSV or TSV) containing updated group
            information.
        merge_column : str
            The column in the metadata file to use for updating the group information.

        Raises
        ------
        ValueError
            If the metadata file is not a valid CSV or TSV file or if the necessary
            columns are missing.
        """
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
        self.sample_metadata["group"] = (
            self.sample_metadata["filename"]
            .map(filename_to_group)
            .fillna(self.sample_metadata["group"])
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
