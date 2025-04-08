# Standard library imports
import os
import pkg_resources
import re
from importlib import resources
from typing import List, Optional, Tuple

# Third-party imports
import pandas as pd


def _load_RDD_metadata(
    external_metadata: Optional[str] = None,
) -> pd.DataFrame:
    """
    Load internal or external reference metadata file containing ontology annotations.

    Parameters
    ----------
    external_metadata : str, optional
        Path to an external metadata file (CSV, TSV, or TXT). If None, internal default metadata is loaded.

    Returns
    -------
    pd.DataFrame
        Reference metadata DataFrame with a cleaned 'filename' column.

    Raises
    ------
    FileNotFoundError
        If the specified external file does not exist.
    ValueError
        If the file format is not supported.
    """
    if external_metadata:
        if not external_metadata.lower().endswith((".csv", ".tsv", ".txt")):
            raise ValueError(
                "External metadata file must be a CSV, TSV, or TXT."
            )
        sep = (
            "\t"
            if external_metadata.lower().endswith((".tsv", ".txt"))
            else ","
        )
        try:
            reference_metadata = pd.read_csv(external_metadata, sep=sep)
        except FileNotFoundError:
            raise FileNotFoundError(
                f"External metadata file '{external_metadata}' not found."
            )
    else:
        try:
            with resources.open_text(
                "data", "foodomics_multiproject_metadata.txt"
            ) as stream:
                reference_metadata = pd.read_csv(stream, sep="\t")
        except (ModuleNotFoundError, ImportError):
            stream = pkg_resources.resource_stream(
                __name__, "data/foodomics_multiproject_metadata.txt"
            )
            reference_metadata = pd.read_csv(stream, sep="\t")

    # process filename for posterior matching
    reference_metadata["filename"] = remove_filename_extension(
        reference_metadata["filename"]
    )
    return reference_metadata


def _load_sample_types(
    reference_metadata: pd.DataFrame,
    simple_complex: str = "all",
    ontology_columns: Optional[List[str]] = None,
) -> Tuple[pd.DataFrame, Optional[List[str]]]:
    """
    Extract and optionally rename ontology columns from the reference metadata.

    Parameters
    ----------
    reference_metadata : pd.DataFrame
        Reference metadata containing ontology and sample_type info.
    simple_complex : str, optional
        Filter for 'simple', 'complex', or 'all' samples (default is 'all').
    ontology_columns : list of str, optional
        Ordered list of ontology columns to be renamed with level suffixes. If None, default ontology structure is used.

    Returns
    -------
    Tuple[pd.DataFrame, list of str or None]
        A tuple containing:
        - A DataFrame with ontology information (indexed by filename).
        - A list of renamed ontology column names, or None if defaults were used.
    """
    if simple_complex != "all":
        reference_metadata = reference_metadata[
            reference_metadata["simple_complex"] == simple_complex
        ]

    if ontology_columns is None:
        sample_type_cols = ["sample_name"] + [
            col
            for col in reference_metadata.columns
            if re.match(r"sample_type_group\d+$", col)
        ]
        df = reference_metadata[["filename", *sample_type_cols]]
        return df.set_index("filename"), None

    # Rename user-defined ontology columns
    renamed_columns = [f"{col}{i+1}" for i, col in enumerate(ontology_columns)]
    renamer = dict(zip(ontology_columns, renamed_columns))
    df = reference_metadata[
        ["filename", "sample_name", *ontology_columns]
    ].rename(columns=renamer)
    return df.set_index("filename"), renamed_columns


def _validate_groups(
    gnps_network: pd.DataFrame, groups_included: List[str]
) -> None:
    """
    Validates that the provided group names exist in the GNPS network data.

    Parameters
    ----------
    gnps_network : pd.DataFrame
        The GNPS network DataFrame to validate against.
    groups_included : List[str]
        The groups to validate.

    Raises
    ------
    ValueError
        If any of the group names in `groups_included` are invalid.
    """
    valid_groups = set(gnps_network["DefaultGroups"].unique())

    # Check included groups
    invalid_included_groups = set(groups_included) - valid_groups
    if len(invalid_included_groups) > 0:
        raise ValueError(
            f"The following groups in groups_included are invalid: {invalid_included_groups}"
        )


def normalize_network(gnps_network, sample_groups=None, reference_groups=None):
    """
    Normalize the GNPS network to extract unique (filename, cluster_index) pairs.

    Parameters
    ----------
    gnps_network : pd.DataFrame
        Raw GNPS network file.
    sample_groups : list of str, optional
        Sample groups to retain.
    reference_groups : list of str, optional
        Reference groups to retain.

    Returns
    -------
    pd.DataFrame
        Normalized DataFrame with columns 'filename' and 'cluster_index'.
    """
    network = gnps_network.copy()

    if "UniqueFileSources" in network.columns:
        groups = {f"G{i}" for i in range(1, 7)}
        groups_excluded = list(
            groups - set([*sample_groups, *reference_groups])
        )
        df_selected = gnps_network[
            (gnps_network[sample_groups] > 0).all(axis=1)
            & (gnps_network[reference_groups] > 0).any(axis=1)
            & (gnps_network[groups_excluded] == 0).all(axis=1)
        ].copy()
        df_exploded = df_selected.assign(
            filename=df_selected["UniqueFileSources"].str.split("|")
        ).explode("filename")

        df_normalized = (
            df_exploded[["filename", "cluster index"]]
            .drop_duplicates()
            .reset_index(drop=True)
        )
        df_normalized.rename(
            columns={"cluster index": "cluster_index"}, inplace=True
        )
        df_normalized["filename"] = remove_filename_extension(
            df_normalized["filename"]
        )
        return df_normalized
    else:
        network["#Filename"] = network["#Filename"].str.replace(
            "input_spectra/", ""
        )
        network.rename(
            columns={"#Filename": "filename", "#ClusterIdx": "cluster_index"},
            inplace=True,
        )
        df_normalized = (
            network[["filename", "cluster_index"]]
            .drop_duplicates()
            .reset_index(drop=True)
        )
        df_normalized["filename"] = remove_filename_extension(
            df_normalized["filename"]
        )

        return df_normalized


def get_sample_metadata(
    raw_gnps_network=None,
    sample_groups=None,
    external_sample_metadata=None,
    filename_col="filename",
):
    """
    Load sample metadata from an external file or extract it from the GNPS network.

    Parameters
    ----------
    raw_gnps_network : pd.DataFrame, optional
        Raw GNPS network DataFrame, required if no external metadata is provided.
    sample_groups : list of str, optional
        List of sample group names to extract from the GNPS network.
    external_sample_metadata : str, optional
        Path to an external sample metadata file (CSV).
    filename_col : str, optional
        Column name for filenames in the metadata file. Default is "filename".

    Returns
    -------
    pd.DataFrame
        Sample metadata DataFrame with columns 'filename' and 'group'.

    Raises
    ------
    ValueError
        If both raw_gnps_network and external_sample_metadata are None.
    """

    if external_sample_metadata:

        if not external_sample_metadata.lower().endswith(
            (".csv", ".tsv", ".txt")
        ):
            raise ValueError(
                "External metadata file must be a CSV, TSV, or TXT."
            )
        sep = (
            "\t"
            if external_sample_metadata.lower().endswith((".tsv", ".txt"))
            else ","
        )
        try:
            sample_metadata = pd.read_csv(external_sample_metadata, sep=sep)
        except FileNotFoundError:
            raise FileNotFoundError(
                f"External metadata file '{external_sample_metadata}' not found."
            )
        sample_metadata.rename(
            columns={filename_col: "filename"}, inplace=True
        )
        sample_metadata["filename"] = remove_filename_extension(
            sample_metadata["filename"]
        )
        return sample_metadata
    else:
        df_filtered = raw_gnps_network[
            ~raw_gnps_network["DefaultGroups"].str.contains(",")
        ]
        df_selected = df_filtered[
            df_filtered["DefaultGroups"].isin(sample_groups)
        ]
        df_exploded_files = df_selected.assign(
            UniqueFileSources=df_selected["UniqueFileSources"].str.split("|")
        ).explode("UniqueFileSources")
        sample_metadata = df_exploded_files[
            ["DefaultGroups", "UniqueFileSources"]
        ].rename(
            columns={"DefaultGroups": "group", "UniqueFileSources": "filename"}
        )
        sample_metadata["filename"] = remove_filename_extension(
            sample_metadata["filename"]
        )
        return sample_metadata.drop_duplicates().reset_index(drop=True)


def split_reference_sample(
    normalized_gnps_network,
    reference_metadata,
    sample_metadata,
    sample_group_col="group",
    reference_name_col="sample_name",
):
    """
    Splits a normalized GNPS network DataFrame into reference and sample subsets
    by matching filenames in `reference_metadata` and `sample_metadata`.

    Parameters
    ----------
    normalized_gnps_network : pd.DataFrame
        A normalized GNPS network DataFrame
    reference_metadata : pd.DataFrame
        A DataFrame containing reference filenames (and other metadata).
        Must include columns like ['filename', reference_name_col].
    sample_metadata : pd.DataFrame
        A DataFrame containing sample filenames (and other metadata).
        Must include columns like ['filename', sample_group_col].
    sample_group_col : str, optional
        The column name in `sample_metadata` indicating sample group,
        by default 'group'.
    reference_name_col : str, optional
        The column name in `reference_metadata` indicating reference name or ID,
        by default 'sample_name'.

    Returns
    -------
    (pd.DataFrame, pd.DataFrame)
        A tuple of (sample_df, reference_df) DataFrames, each a subset of
        the noramlized gnps network. The first contains rows whose filenames match
        those in the sample metadata, and the second contains rows whose filenames
        match those in the reference metadata.
    """

    sample_clusters = (
        pd.merge(
            normalized_gnps_network,
            sample_metadata[["filename", sample_group_col]],
            on="filename",
            how="inner",
        )
        .drop_duplicates()
        .reset_index(drop=True)
    )

    reference_clusters = (
        pd.merge(
            normalized_gnps_network,
            reference_metadata[["filename", reference_name_col]],
            on="filename",
            how="inner",
        )
        .drop_duplicates()
        .reset_index(drop=True)
    )

    return sample_clusters, reference_clusters


def remove_filename_extension(filename_col):
    """
    Removes the file extension from a filename column in a DataFrame.

    Parameters
    ----------
    filename_col : str
        The name of the column containing filenames.

    Returns
    -------
    str
        The filename without its extension.
    """
    return filename_col.str.replace(r"\.[^.]+$", "", regex=True)


def RDD_counts_to_wide(
    RDD_counts: pd.DataFrame, level: int = None
) -> pd.DataFrame:
    """
    Convert the RDD counts dataframe from long to wide format for a specific
    ontology level, with 'group' as part of the columns. If the data is already
    filtered by level, no level needs to be passed.

    Parameters
    ----------
    RDD_counts : pd.DataFrame
        A long-format DataFrame with columns ['filename', 'reference_type', 'count',
        'level', 'group'] representing RDD counts across different ontology
        levels and groups.
    level : int, optional
        The ontology level to filter by before converting to wide format. If
        None, the level is inferred from the data. Defaults to None.

    Returns
    -------
    pd.DataFrame
        A wide-format DataFrame where each combination of 'reference_type' and
        'group' becomes a column, and rows are indexed by 'filename'.

    Raises
    ------
    ValueError
        If multiple levels are found in the data and `level` is not specified,
        or if no data is available for the specified level.
    """
    # If level is not specified, infer the level from the data
    if level is None:
        levels_in_data = RDD_counts["level"].unique()
        if len(levels_in_data) > 1:
            raise ValueError(
                "Multiple levels found in the data. Please specify a level to convert to wide format."
            )
        level = levels_in_data[
            0
        ]  # If the data is already filtered, use that level

    # Filter the RDD counts dataframe by the specified level, if not already filtered
    filtered_RDD_counts = RDD_counts[RDD_counts["level"] == level]

    if filtered_RDD_counts.empty:
        raise ValueError(f"No data available for level {level}")

    # Pivot the filtered dataframe to wide format with 'reference_type' and 'group' as columns
    RDD_counts_wide = filtered_RDD_counts.pivot_table(
        index="filename",
        columns="reference_type",
        values="count",
        fill_value=0,
    )

    # Convert only numerical columns (count values) to int
    RDD_counts_wide = RDD_counts_wide.apply(
        lambda col: col.astype(int) if col.dtype == "float64" else col
    )
    group_df = (
        filtered_RDD_counts[["filename", "group"]]
        .drop_duplicates()
        .set_index("filename")
    )
    wide_format_RDD_counts = RDD_counts_wide.join(group_df)

    return wide_format_RDD_counts


def calculate_proportions(
    RDD_counts: pd.DataFrame, level: int = None
) -> pd.DataFrame:
    """
    Calculate the proportion of each reference type within each sample for a given
    level.

    Parameters
    ----------
    RDD_counts : pd.DataFrame
        A long-format DataFrame with columns ['filename', 'reference_type', 'count',
        'level', 'group'] representing RDD counts across different ontology
        levels and groups.
    level : int, optional
        The ontology level to filter by before calculating proportions. If None,
        the level is inferred from the data. Defaults to None.

    Returns
    -------
    pd.DataFrame
        A wide-format DataFrame where each reference type column contains the
        proportion of that reference type within each sample (row). Rows are indexed
        by 'filename', and proportions sum to 1 for each sample.

    Raises
    ------
    ValueError
        If multiple levels are found in the data and `level` is not specified.
    """
    # If level is not specified, infer the level from the data
    if level is None:
        levels_in_data = RDD_counts["level"].unique()
        if len(levels_in_data) > 1:
            raise ValueError(
                "Multiple levels found in the data. Please specify a level to calculate proportions."
            )
        level = levels_in_data[
            0
        ]  # If the data is already filtered, use that level

    # Use the existing function to convert to wide format
    df_wide = RDD_counts_to_wide(RDD_counts, level)

    # Identify numeric columns (reference type counts)
    numeric_cols = df_wide.select_dtypes(include=[float, int]).columns

    # Calculate proportions across reference types (columns) for each sample (row)
    df_proportions = df_wide.copy()
    df_proportions[numeric_cols] = (
        df_proportions[numeric_cols]
        .div(df_proportions[numeric_cols].sum(axis=1), axis=0)
        .fillna(0)
    )

    return df_proportions
