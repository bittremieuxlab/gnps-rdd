import os
import pytest
import pandas as pd
import sys
# Add the path to the 'rdd' directory to the system path
project_path = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "rdd")
)
sys.path.append(project_path)

from utils import (
    _load_RDD_metadata,
    _load_sample_types,
    _validate_groups,
    normalize_network,
    get_sample_metadata,
    split_reference_sample,
    remove_filename_extension,
    RDD_counts_to_wide,
    calculate_proportions
)


@pytest.fixture
def sample_reference_metadata():
    return pd.DataFrame({
        "filename": ["ref1.mzXML", "ref2.mzXML"],
        "simple_complex": ["simple", "complex"],
        "sample_name": ["apple", "meat"],
        "sample_type_group1": ["fruit", "animal"],
        "sample_type_group2": ["fruit", "meat"]
    })


@pytest.fixture
def sample_gnps_network():
    return pd.DataFrame({
        "DefaultGroups": ["G1", "G2"],
        "UniqueFileSources": ["file1.mzXML|file2.mzXML", "file3.mzXML"],
        "G1": [1, 0],
        "G2": [0, 1],
        "G3": [0, 0],
        "G4": [0, 0],
        "G5": [0, 0],
        "G6": [0, 0],
        "cluster index": [101, 102]
    })


def test_load_sample_types_default_structure(sample_reference_metadata):
    result, renamed = _load_sample_types(sample_reference_metadata)
    assert result.index.name == "filename"
    assert renamed is None
    assert "sample_type_group1" in result.columns


def test_load_sample_types_custom_columns(sample_reference_metadata):
    result, renamed = _load_sample_types(
        sample_reference_metadata, ontology_columns=["sample_name", "sample_type_group1"]
    )
    assert "sample_name1" in result.columns
    assert "sample_type_group12" in result.columns
    assert renamed == ["sample_name1", "sample_type_group12"]


def test_validate_groups_pass(sample_gnps_network):
    _validate_groups(sample_gnps_network, ["G1", "G2"])  # Should not raise


def test_validate_groups_fail(sample_gnps_network):
    with pytest.raises(ValueError):
        _validate_groups(sample_gnps_network, ["G3"])


def test_remove_filename_extension():
    filenames = pd.Series(["sample1.mzXML", "sample2.mzML"])
    result = remove_filename_extension(filenames)
    assert result.tolist() == ["sample1", "sample2"]


def test_normalize_network_network(sample_gnps_network):
    df = normalize_network(sample_gnps_network, sample_groups=["G1"], reference_groups=["G2"])
    assert {"filename", "cluster_index"}.issubset(df.columns)


def test_split_reference_sample():
    normalized = pd.DataFrame({
        "filename": ["sample1", "ref1"],
        "cluster_index": [1, 2]
    })

    sample_meta = pd.DataFrame({
        "filename": ["sample1"],
        "group": ["G1"]
    })

    ref_meta = pd.DataFrame({
        "filename": ["ref1"],
        "sample_name": ["apple"]
    })

    sample_clusters, ref_clusters = split_reference_sample(normalized, ref_meta, sample_meta)
    assert not sample_clusters.empty
    assert not ref_clusters.empty
    assert "cluster_index" in sample_clusters.columns


def test_RDD_counts_to_wide_and_proportions():
    df = pd.DataFrame({
        "filename": ["f1", "f2", "f1", "f2"],
        "reference_type": ["r1", "r1", "r2", "r2"],
        "count": [10, 20, 30, 40],
        "level": [1, 1, 1, 1],
        "group": ["G1", "G1", "G1", "G1"]
    })

    wide = RDD_counts_to_wide(df, level=1)
    assert "r1" in wide.columns and "r2" in wide.columns
    assert wide.index.name == "filename"

    proportions = calculate_proportions(df, level=1)
    assert abs(proportions.loc["f1", "r1"] + proportions.loc["f1", "r2"] - 1) < 1e-6
