import os
import pytest
import pandas as pd
import sys

# Add the path to the repository root to the system path
project_path = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if project_path not in sys.path:
    sys.path.insert(0, project_path)

from rdd.utils import (
    _load_RDD_metadata,
    _load_sample_types,
    _validate_groups,
    normalize_network,
    get_sample_metadata,
    split_reference_sample,
    remove_filename_extension,
    RDD_counts_to_wide,
    calculate_proportions,
)


@pytest.fixture
def sample_reference_metadata():
    return pd.DataFrame(
        {
            "filename": ["ref1.mzXML", "ref2.mzXML"],
            "simple_complex": ["simple", "complex"],
            "sample_name": ["apple", "meat"],
            "sample_type_group1": ["fruit", "animal"],
            "sample_type_group2": ["fruit", "meat"],
        }
    )


@pytest.fixture
def sample_gnps_network():
    return pd.DataFrame(
        {
            "DefaultGroups": ["G1", "G2"],
            "UniqueFileSources": ["file1.mzXML|file2.mzXML", "file3.mzXML"],
            "G1": [1, 0],
            "G2": [0, 1],
            "G3": [0, 0],
            "G4": [0, 0],
            "G5": [0, 0],
            "G6": [0, 0],
            "cluster index": [101, 102],
        }
    )


def test_load_sample_types_default_structure(sample_reference_metadata):
    result, _ = _load_sample_types(sample_reference_metadata)
    assert result.index.name == "filename"
    assert "sample_type_group1" in result.columns


def test_load_sample_types_custom_columns(sample_reference_metadata):
    result, renamed = _load_sample_types(
        sample_reference_metadata,
        ontology_columns=["sample_name", "sample_type_group1"],
    )
    assert "sample_name1" in result.columns
    assert "sample_type_group12" in result.columns
    assert renamed == ["sample_name1", "sample_type_group12"]


def test_load_sample_types_with_simple_filter(sample_reference_metadata):
    """Test _load_sample_types with simple_complex filter."""
    result, _ = _load_sample_types(
        sample_reference_metadata,
        simple_complex="simple",
    )
    # Should only have simple samples
    assert len(result) == 1  # Only ref1 is simple


def test_load_sample_types_with_complex_filter(sample_reference_metadata):
    """Test _load_sample_types with simple_complex filter for complex."""
    result, _ = _load_sample_types(
        sample_reference_metadata,
        simple_complex="complex",
    )
    # Should only have complex samples
    assert len(result) == 1  # Only ref2 is complex


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
    df = normalize_network(
        sample_gnps_network, sample_groups=["G1"], reference_groups=["G2"]
    )
    assert {"filename", "cluster_index"}.issubset(df.columns)


def test_split_reference_sample():
    normalized = pd.DataFrame(
        {"filename": ["sample1", "ref1"], "cluster_index": [1, 2]}
    )

    sample_meta = pd.DataFrame({"filename": ["sample1"], "group": ["G1"]})

    ref_meta = pd.DataFrame({"filename": ["ref1"], "sample_name": ["apple"]})

    sample_clusters, ref_clusters = split_reference_sample(
        normalized, ref_meta, sample_meta
    )
    assert not sample_clusters.empty
    assert not ref_clusters.empty
    assert "cluster_index" in sample_clusters.columns


def test_RDD_counts_to_wide_and_proportions():
    df = pd.DataFrame(
        {
            "filename": ["f1", "f2", "f1", "f2"],
            "reference_type": ["r1", "r1", "r2", "r2"],
            "count": [10, 20, 30, 40],
            "level": [1, 1, 1, 1],
            "group": ["G1", "G1", "G1", "G1"],
        }
    )

    wide = RDD_counts_to_wide(df, level=1)
    assert "r1" in wide.columns and "r2" in wide.columns
    assert wide.index.name == "filename"

    proportions = calculate_proportions(df, level=1)
    assert (
        abs(proportions.loc["f1", "r1"] + proportions.loc["f1", "r2"] - 1)
        < 1e-6
    )


# ============================================================================
# Additional get_sample_metadata Tests
# ============================================================================


def test_get_sample_metadata_without_inputs():
    """Test get_sample_metadata raises error when neither input is provided."""
    with pytest.raises(ValueError, match="raw_gnps_network is required"):
        get_sample_metadata(
            raw_gnps_network=None,
            sample_groups=None,
            external_sample_metadata=None,
        )


def test_get_sample_metadata_without_sample_groups():
    """Test get_sample_metadata raises error when sample_groups not provided."""
    gnps_data = pd.DataFrame(
        {
            "cluster index": [1],
            "DefaultGroups": ["G1"],
            "UniqueFileSources": ["file1.mzML"],
        }
    )

    with pytest.raises(ValueError, match="sample_groups must be provided"):
        get_sample_metadata(
            raw_gnps_network=gnps_data,
            sample_groups=None,
            external_sample_metadata=None,
        )


def test_get_sample_metadata_with_invalid_file_format(tmp_path):
    """Test get_sample_metadata with invalid file format."""
    invalid_file = tmp_path / "metadata.json"
    invalid_file.write_text('{"test": "data"}')

    with pytest.raises(
        ValueError, match="External metadata file must be a CSV, TSV, or TXT"
    ):
        get_sample_metadata(
            external_sample_metadata=str(invalid_file),
        )


def test_get_sample_metadata_with_missing_file():
    """Test get_sample_metadata with non-existent file."""
    with pytest.raises(FileNotFoundError, match="not found"):
        get_sample_metadata(
            external_sample_metadata="/nonexistent/path/metadata.csv",
        )


def test_get_sample_metadata_with_invalid_column(tmp_path):
    """Test get_sample_metadata with missing filename column."""
    metadata = """sample_id,group
sample1,G1
sample2,G2"""

    metadata_file = tmp_path / "metadata.csv"
    metadata_file.write_text(metadata)

    with pytest.raises(KeyError, match="Column 'nonexistent_col' not found"):
        get_sample_metadata(
            external_sample_metadata=str(metadata_file),
            filename_col="nonexistent_col",
        )


def test_get_sample_metadata_with_tsv_file(tmp_path):
    """Test get_sample_metadata correctly handles TSV files."""
    metadata = "filename\tgroup\nfile1\tG1\nfile2\tG2"

    metadata_file = tmp_path / "metadata.tsv"
    metadata_file.write_text(metadata)

    result = get_sample_metadata(external_sample_metadata=str(metadata_file))

    assert len(result) == 2
    assert "filename" in result.columns
    assert "group" in result.columns


def test_get_sample_metadata_with_txt_file(tmp_path):
    """Test get_sample_metadata correctly handles TXT files."""
    metadata = "filename\tgroup\nfile1\tG1\nfile2\tG2"

    metadata_file = tmp_path / "metadata.txt"
    metadata_file.write_text(metadata)

    result = get_sample_metadata(external_sample_metadata=str(metadata_file))

    assert len(result) == 2
    assert "filename" in result.columns
    assert "group" in result.columns


# ============================================================================
# Additional _load_RDD_metadata Tests
# ============================================================================


def test_load_rdd_metadata_with_invalid_format(tmp_path):
    """Test _load_RDD_metadata with invalid file format."""
    invalid_file = tmp_path / "metadata.xml"
    invalid_file.write_text("<xml>test</xml>")

    with pytest.raises(
        ValueError, match="External metadata file must be a CSV, TSV, or TXT"
    ):
        _load_RDD_metadata(external_metadata=str(invalid_file))


def test_load_rdd_metadata_with_missing_file():
    """Test _load_RDD_metadata with non-existent file."""
    with pytest.raises(FileNotFoundError, match="not found"):
        _load_RDD_metadata(external_metadata="/nonexistent/path/metadata.csv")


def test_load_rdd_metadata_with_tsv_file(tmp_path):
    """Test _load_RDD_metadata correctly handles TSV files."""
    metadata = (
        "filename\tsample_name\tsample_type_group1\nfile1.mzML\tref1\ttype1"
    )

    metadata_file = tmp_path / "metadata.tsv"
    metadata_file.write_text(metadata)

    result = _load_RDD_metadata(external_metadata=str(metadata_file))

    assert len(result) == 1
    assert "filename" in result.columns
    # Filename extension should be removed
    assert result["filename"].iloc[0] == "file1"


# ============================================================================
# Additional calculate_proportions Tests
# ============================================================================


def test_calculate_proportions_with_multiple_levels():
    """Test calculate_proportions raises error with multiple levels."""
    rdd_counts = pd.DataFrame(
        {
            "filename": ["file1", "file1", "file2", "file2"],
            "reference_type": ["type1", "type2", "type1", "type2"],
            "count": [10, 20, 15, 25],
            "level": [1, 1, 2, 2],  # Multiple levels
            "group": ["G1", "G1", "G2", "G2"],
        }
    )

    with pytest.raises(ValueError, match="Multiple levels found"):
        calculate_proportions(rdd_counts)


def test_calculate_proportions_with_explicit_level():
    """Test calculate_proportions with explicit level parameter."""
    rdd_counts = pd.DataFrame(
        {
            "filename": ["file1", "file1", "file2", "file2"],
            "reference_type": ["type1", "type2", "type1", "type2"],
            "count": [10, 20, 15, 25],
            "level": [1, 1, 2, 2],
            "group": ["G1", "G1", "G2", "G2"],
        }
    )

    # Should work with explicit level
    result = calculate_proportions(rdd_counts, level=1)

    assert "type1" in result.columns
    assert "type2" in result.columns
    assert result.loc["file1", "type1"] == pytest.approx(10 / 30)
    assert result.loc["file1", "type2"] == pytest.approx(20 / 30)
