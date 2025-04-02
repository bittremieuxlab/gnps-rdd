# RDDcounts.py

# Standard library imports
import os
from typing import List, Optional, Tuple

# Third-party imports
import numpy as np
import pandas as pd

# Internal imports
from utils import _load_RDD_metadata, _load_sample_types, _validate_groups, get_sample_metadata, normalize_network, split_reference_sample

class RDDCounts:
    def __init__(
        self,
        gnps_network_path: str,
        sample_types: str,
        sample_groups: List[str] = None,
        reference_groups: List[str] = None,
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

        self.reference_metadata = _load_RDD_metadata(external_reference_metadata)
        self.sample_metadata = get_sample_metadata(
            self.raw_gnps_network,
            sample_groups,
            external_sample_metadata,
            filename_col="filename"
        )
        if ontology_columns is not None:
            self.ontology_columns_renamed = [f"{col}{i+1}" for i, col in enumerate(ontology_columns)]
            self.sample_types_df = _load_sample_types(
                self.reference_metadata,
                sample_types,
                ontology_columns=ontology_columns
            )
        else:
            self.ontology_columns_renamed = None
            self.sample_types_df = _load_sample_types(
                self.reference_metadata,
                sample_types
            )
        self.normalized_network = normalize_network(self.raw_gnps_network,self.sample_groups, self.reference_groups)
        self.file_level_counts = self.file_counts(self.normalized_network, self.reference_metadata, self.sample_metadata)
        self.counts = self.create_RDD_counts_all_levels()
       

    def file_counts(self, normalized_gnps_network, reference_metadata, sample_metadata, sample_group_col="group", reference_name_col="sample_name"):
        """
        Returns
        -------
        pd.DataFrame
            A DataFrame containing file-level counts.
        """
        sample_clusters, reference_clusters = split_reference_sample(normalized_gnps_network, reference_metadata, sample_metadata, sample_group_col, reference_name_col)        
        shared_clusters = reference_clusters.merge(sample_clusters, on="cluster_index", suffixes=("_reference", "_sample"))
        cluster_count = shared_clusters.groupby(["filename_sample", reference_name_col]).size().unstack(fill_value=0)
        cluster_count = cluster_count.drop_duplicates()
        cluster_count_long = cluster_count.stack().reset_index(name="count")
        cluster_count_long.rename(columns={reference_name_col: "reference_type", "filename_sample": "filename"}, inplace=True)
        cluster_count_long["level"] = 0
        cluster_count_long["group"] = cluster_count_long.merge(self.sample_metadata, on="filename", how="inner")["group"]
        return cluster_count_long

    def create_RDD_counts_all_levels(self) -> pd.DataFrame:
        """
        Generates RDD counts across all ontology levels and compiles them into a
        single DataFrame. This function creates level-specific counts by grouping the
        file-level counts and applies a filter for samples appearing less frequent than
        water.Returns the data in long format.

        Returns
        -------
        pd.DataFrame
            A concatenated DataFrame containing RDD counts across all levels,
            including 'filename', 'reference_type', 'count', 'level', and 'group' columns.
        """
        RDD_counts_file_level = self.file_level_counts
        RDD_counts_all_levels = [
            RDD_counts_file_level
        ]  # Initialize a list for storing data at all levels
        RDD_counts_file_level_sample_types = RDD_counts_file_level.merge(
            self.sample_types_df, left_on="reference_type", right_on="sample_name"
        ).drop_duplicates()

        sample_metadata_map = self.sample_metadata.set_index("filename")[
            "group"
        ].to_dict()  # Create the mapping once

        for level in range(1, self.levels + 1):
            if self.ontology_columns_renamed:
                ontology_col = self.ontology_columns_renamed[level - 1]
            else:
                ontology_col = f"sample_type_group{level}"

            RDD_counts_level = (
                RDD_counts_file_level_sample_types.groupby(["filename", ontology_col])["count"]
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
                columns_to_modify = wide_format_counts.columns.difference(["filename", "water"])
                wide_format_counts.loc[:, columns_to_modify] = wide_format_counts.loc[:, columns_to_modify].where(
                    wide_format_counts.loc[:, columns_to_modify].gt(water_counts, axis=0),
                    0,
                )
                wide_format_counts = wide_format_counts.drop(columns=["water"])

            wide_format_counts = wide_format_counts.loc[:, (wide_format_counts != 0).any(axis=0)]

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
            # Set source and target level for the current iteration
            source_level = i
            target_level = i + 1

            # Group RDD counts at target level
            target_counts = (
                counts[counts["level"] == target_level]
                .groupby("reference_type")["count"]
                .sum()
                .reset_index()
            )

            # Merge with sample_types to find corresponding source level
            merged_df = pd.merge(
                target_counts,
                self.sample_types,
                left_on="reference_type",
                right_on=f"sample_type_group{target_level}",
            )[
                [f"sample_type_group{source_level}", "reference_type", "count"]
            ].drop_duplicates()

            # Rename columns for source-target relationship and add levels for uniqueness
            flow = merged_df.rename(
                columns={
                    f"sample_type_group{source_level}": "source",
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
        all_nodes = pd.concat(
            [flows_df["source"], flows_df["target"]]
        ).unique()
        processes_df = pd.DataFrame(
            {
                "id": all_nodes,
                "level": [int(node.split("_")[-1]) for node in all_nodes],
            }
        ).set_index("id")

        return flows_df, processes_df
