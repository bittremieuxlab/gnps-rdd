import os
import sys
import pytest
import pandas as pd

# Add the path to the repository root to the system path
project_path = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if project_path not in sys.path:
    sys.path.insert(0, project_path)

from rdd.RDDcounts import RDDCounts


@pytest.mark.parametrize("fixture_type", ["gnps1", "gnps2"])
def test_rddcounts_initialization(load_test_files, fixture_type):
    gnps_path, sample_meta, ref_meta = load_test_files(fixture_type)

    # Condition for direct vs. network initialization
    if fixture_type == "gnps1":
        sample_groups = ["G1"]
        reference_groups = ["G4"]
        sample_meta = None
    else:
        sample_groups = None  # No sample groups for direct format
        reference_groups = None  # No reference groups for direct format
        sample_meta = str(sample_meta)

    rdd = RDDCounts(
        gnps_network_path=str(gnps_path),
        sample_types="all",
        sample_groups=sample_groups,
        reference_groups=reference_groups,
        levels=6,
        external_reference_metadata=str(ref_meta),
        external_sample_metadata=sample_meta,
    )

    assert isinstance(rdd.raw_gnps_network, pd.DataFrame)
    assert isinstance(rdd.counts, pd.DataFrame)
    assert not rdd.counts.empty


@pytest.mark.parametrize("fixture_type", ["gnps1", "gnps2"])
def test_filter_counts_functionality(load_test_files, fixture_type):
    gnps_path, sample_meta, ref_meta = load_test_files(fixture_type)
    rdd = RDDCounts(
        gnps_network_path=str(gnps_path),
        sample_types="all",
        sample_groups=["G1"],
        reference_groups=["G4"],
        levels=3,
        external_reference_metadata=str(ref_meta),
        external_sample_metadata=str(sample_meta),
    )

    filtered = rdd.filter_counts(level=1)
    assert "reference_type" in filtered.columns
    assert filtered["level"].eq(1).all()


@pytest.mark.parametrize("fixture_type", ["gnps1", "gnps2"])
def test_generate_flows_structure(load_test_files, fixture_type):
    gnps_path, sample_meta, ref_meta = load_test_files(fixture_type)
    rdd = RDDCounts(
        gnps_network_path=str(gnps_path),
        sample_types="all",
        sample_groups=["G1"],
        reference_groups=["G4"],
        levels=3,
        external_reference_metadata=str(ref_meta),
        external_sample_metadata=str(sample_meta),
    )

    flows_df, processes_df = rdd.generate_RDDflows()
    assert {"source", "target", "value"}.issubset(flows_df.columns)
    assert "level" in processes_df.columns


# ============================================================================
# Validation and Error Handling Tests
# ============================================================================


def test_rddcounts_with_task_id_validation():
    """Test that providing both task_id and gnps_network_path raises error."""
    with pytest.raises(ValueError, match="Provide exactly one"):
        RDDCounts(
            sample_types="all",
            task_id="fake_task_id",
            gnps_network_path="fake_path.tsv",
            sample_groups=["G1"],
            reference_groups=["G2"],
        )


def test_rddcounts_with_neither_task_nor_path():
    """Test that providing neither task_id nor gnps_network_path raises error."""
    with pytest.raises(ValueError, match="Provide exactly one"):
        RDDCounts(
            sample_types="all",
            task_id=None,
            gnps_network_path=None,
            sample_groups=["G1"],
            reference_groups=["G2"],
        )


def test_rddcounts_levels_exceeds_ontology_columns(tmp_path):
    """Test that levels exceeding ontology columns raises error."""
    # Create minimal test data
    gnps_data = """cluster index\tsum(precursor intensity)\tRTMean\tG1\tG2\tDefaultGroups\tUniqueFileSources
1\t1000.0\t1.5\t5\t0\tG1\tfile1.mzML|file2.mzML
2\t2000.0\t2.5\t0\t3\tG2\tfile3.mzML"""

    gnps_file = tmp_path / "gnps_test.tsv"
    gnps_file.write_text(gnps_data)

    # Create metadata with only 2 ontology columns
    ref_metadata = """filename,sample_name,ontology1,ontology2,simple_complex
file1,ref1,type1,subtype1,simple
file2,ref2,type2,subtype2,simple"""

    ref_file = tmp_path / "ref_metadata.csv"
    ref_file.write_text(ref_metadata)

    with pytest.raises(
        ValueError, match="levels .* exceeds provided ontology columns"
    ):
        RDDCounts(
            sample_types="all",
            gnps_network_path=str(gnps_file),
            sample_groups=["G1"],
            reference_groups=["G2"],
            levels=5,  # More than 2 ontology columns
            ontology_columns=["ontology1", "ontology2"],
            external_reference_metadata=str(ref_file),
        )


def test_get_ontology_column_for_level_invalid_range(load_test_files):
    """Test that invalid level raises error in _get_ontology_column_for_level."""
    gnps_data, sample_meta, ref_meta = load_test_files("gnps1")

    rdd = RDDCounts(
        sample_types="all",
        gnps_network_path=str(gnps_data),
        sample_groups=["G1", "G2"],
        reference_groups=["G3"],
        external_reference_metadata=str(ref_meta),
        external_sample_metadata=str(sample_meta),
        ontology_columns=["sample_type_group1", "sample_type_group2"],
        levels=2,
    )

    # Test level too low
    with pytest.raises(ValueError, match="Invalid level 0"):
        rdd._get_ontology_column_for_level(0)

    # Test level too high
    with pytest.raises(ValueError, match="Invalid level 3"):
        rdd._get_ontology_column_for_level(3)


def test_rddcounts_get_ontology_column_default_format(load_test_files):
    """Test _get_ontology_column_for_level with default column format."""
    gnps_data, sample_meta, ref_meta = load_test_files("gnps1")

    rdd = RDDCounts(
        sample_types="all",
        gnps_network_path=str(gnps_data),
        sample_groups=["G1", "G2"],
        reference_groups=["G3"],
        external_reference_metadata=str(ref_meta),
        external_sample_metadata=str(sample_meta),
        levels=3,
    )

    # Should return default format when ontology_columns not provided
    assert rdd._get_ontology_column_for_level(1) == "sample_type_group1"
    assert rdd._get_ontology_column_for_level(2) == "sample_type_group2"
    assert rdd._get_ontology_column_for_level(3) == "sample_type_group3"


# ============================================================================
# Filter Counts Tests with Edge Cases
# ============================================================================


def test_filter_counts_without_initialization():
    """Test filter_counts raises error when counts not initialized."""
    # Create a minimal RDDCounts instance
    gnps_data = pd.DataFrame(
        {
            "cluster index": [1],
            "#Filename": ["file1.mzML"],
            "#ClusterIdx": [1],
        }
    )

    # We need to manually create an instance without calling methods that initialize counts
    # This is tricky, so we'll use a mock or set counts to None
    with pytest.raises(ValueError, match="RDD counts are not initialized"):
        rdd = object.__new__(RDDCounts)
        rdd.counts = None
        rdd.filter_counts(level=1)


def test_filter_counts_with_upper_lower_level_validation(load_test_files):
    """Test filter_counts validates upper_level < lower_level."""
    gnps_data, sample_meta, ref_meta = load_test_files("gnps1")

    rdd = RDDCounts(
        sample_types="all",
        gnps_network_path=str(gnps_data),
        sample_groups=["G1", "G2"],
        reference_groups=["G3"],
        external_reference_metadata=str(ref_meta),
        external_sample_metadata=str(sample_meta),
        levels=3,
    )

    # upper_level must be less than lower_level
    with pytest.raises(
        ValueError, match="upper_level must be lower than lower_level"
    ):
        rdd.filter_counts(level=2, upper_level=2, lower_level=1)


def test_filter_counts_requires_upper_level_reference_types(load_test_files):
    """Test filter_counts requires upper_level_reference_types when using upper_level."""
    gnps_data, sample_meta, ref_meta = load_test_files("gnps1")

    rdd = RDDCounts(
        sample_types="all",
        gnps_network_path=str(gnps_data),
        sample_groups=["G1", "G2"],
        reference_groups=["G3"],
        external_reference_metadata=str(ref_meta),
        external_sample_metadata=str(sample_meta),
        levels=3,
    )

    with pytest.raises(
        ValueError, match="Must provide upper_level_reference_types"
    ):
        rdd.filter_counts(
            level=2,
            upper_level=1,
            lower_level=2,
            upper_level_reference_types=None,
        )


def test_filter_counts_with_upper_lower_level_filtering(load_test_files):
    """Test filter_counts with upper and lower level filtering."""
    gnps_data, sample_meta, ref_meta = load_test_files("gnps1")

    rdd = RDDCounts(
        sample_types="all",
        gnps_network_path=str(gnps_data),
        sample_groups=["G1", "G2"],
        reference_groups=["G3"],
        external_reference_metadata=str(ref_meta),
        external_sample_metadata=str(sample_meta),
        levels=3,
    )

    # Get reference types from level 1 first
    level_1_types = (
        rdd.filter_counts(level=1)["reference_type"].unique().tolist()
    )

    if len(level_1_types) > 0:
        # Filter with upper/lower level
        filtered = rdd.filter_counts(
            level=2,
            upper_level=1,
            lower_level=2,
            upper_level_reference_types=level_1_types[:1],  # Use first type
        )

        # Should be filtered to level 2
        assert filtered["level"].eq(2).all()


def test_filter_counts_with_sample_names(load_test_files):
    """Test filter_counts with sample_names filtering."""
    gnps_data, sample_meta, ref_meta = load_test_files("gnps1")

    rdd = RDDCounts(
        sample_types="all",
        gnps_network_path=str(gnps_data),
        sample_groups=["G1", "G2"],
        reference_groups=["G3"],
        external_reference_metadata=str(ref_meta),
        external_sample_metadata=str(sample_meta),
        levels=3,
    )

    # Get available filenames
    all_files = rdd.counts["filename"].unique().tolist()

    if len(all_files) > 0:
        # Filter by specific sample
        filtered = rdd.filter_counts(level=1, sample_names=[all_files[0]])

        # Should only contain the specified sample
        assert filtered["filename"].unique().tolist() == [all_files[0]]


def test_filter_counts_with_reference_types(load_test_files):
    """Test filter_counts with reference_types filtering."""
    gnps_data, sample_meta, ref_meta = load_test_files("gnps1")

    rdd = RDDCounts(
        sample_types="all",
        gnps_network_path=str(gnps_data),
        sample_groups=["G1", "G2"],
        reference_groups=["G3"],
        external_reference_metadata=str(ref_meta),
        external_sample_metadata=str(sample_meta),
        levels=3,
    )

    # Get available reference types
    all_types = rdd.counts["reference_type"].unique().tolist()

    if len(all_types) > 0:
        # Filter by specific reference type
        filtered = rdd.filter_counts(level=1, reference_types=[all_types[0]])

        # Should only contain the specified reference type
        assert filtered["reference_type"].unique().tolist() == [all_types[0]]


def test_filter_counts_with_top_and_reference_types(load_test_files):
    """Test filter_counts with both top_n and explicit reference_types."""
    gnps_path, sample_meta, ref_meta = load_test_files("gnps1")

    rdd = RDDCounts(
        gnps_network_path=str(gnps_path),
        sample_types="all",
        sample_groups=["G1"],
        reference_groups=["G4"],
        external_reference_metadata=str(ref_meta),
        external_sample_metadata=str(sample_meta),
        levels=3,
    )

    # Get available reference types
    all_types = rdd.counts["reference_type"].unique().tolist()

    if len(all_types) > 1:
        # Filter with top_n and explicit reference_types
        filtered = rdd.filter_counts(
            level=1, top_n=10, reference_types=[all_types[0]]
        )

        # Should only contain the specified reference type
        assert filtered["reference_type"].unique().tolist() == [all_types[0]]


def test_generate_flows_with_filename_filter(load_test_files):
    """Test generate_RDDflows with filename_filter parameter."""
    gnps_path, sample_meta, ref_meta = load_test_files("gnps1")

    rdd = RDDCounts(
        gnps_network_path=str(gnps_path),
        sample_types="all",
        sample_groups=["G1"],
        reference_groups=["G4"],
        external_reference_metadata=str(ref_meta),
        external_sample_metadata=str(sample_meta),
        levels=3,
    )

    # Get a valid filename from the counts
    if not rdd.counts.empty:
        test_filename = rdd.counts["filename"].iloc[0]
        flows_df, processes_df = rdd.generate_RDDflows(
            filename_filter=test_filename
        )

        assert isinstance(flows_df, pd.DataFrame)
        assert isinstance(processes_df, pd.DataFrame)


def test_generate_flows_without_renamed_columns(load_test_files):
    """Test generate_RDDflows when ontology_columns_renamed is None."""
    gnps_path, sample_meta, ref_meta = load_test_files("gnps2")

    rdd = RDDCounts(
        gnps_network_path=str(gnps_path),
        sample_types="all",
        external_reference_metadata=str(ref_meta),
        external_sample_metadata=str(sample_meta),
        levels=3,
    )

    # Ensure ontology_columns_renamed is not set
    rdd.ontology_columns_renamed = None

    flows_df, processes_df = rdd.generate_RDDflows()

    assert isinstance(flows_df, pd.DataFrame)
    assert isinstance(processes_df, pd.DataFrame)
    assert {"source", "target", "value"}.issubset(flows_df.columns)
