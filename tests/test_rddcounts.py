import os
import sys
import pytest
import pandas as pd
# Add the path to the 'rdd' directory to the system path
project_path = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "rdd")
)
sys.path.append(project_path)

from rdd.RDDcounts import RDDCounts

@pytest.mark.parametrize("fixture_type", ["gnps1", "gnps2"])
def test_rddcounts_initialization(load_test_files, fixture_type):
    gnps_path, sample_meta, ref_meta = load_test_files(fixture_type)

    # Condition for direct vs. network initialization
    if fixture_type == "gnps1":
        sample_groups = ["G1"]
        reference_groups = ["G4"]
        sample_meta=None
    else:
        sample_groups = None  # No sample groups for direct format
        reference_groups = None  # No reference groups for direct format
        sample_meta= str(sample_meta)

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