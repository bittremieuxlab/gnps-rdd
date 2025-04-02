# conftest.py
import pytest
import pandas as pd
from pathlib import Path


@pytest.fixture
def mock_network_based_gnps(tmp_path):
    """Fixture for the GNPS file with 'DefaultGroups' and 'UniqueFileSources'."""
    # Create mock GNPS file
    data = {
        "DefaultGroups": ["G1", "G4", "G1,G4", "G1,G4", "G1,G4", "G1,G4"],
        "UniqueFileSources": [
            "file_samp_1.mzXML|file_samp_2.mzXML",
            "file_ref1.mzXML|file_ref2.mzXML",
            "file_samp_1.mzXML|file_ref1.mzXML",
            "file_samp_2.mzXML|file_ref2.mzXML",
            "file_samp_1.mzXML|file_ref3.mzXML",
            "file_samp_2.mzXML|file_ref3.mzXML",
        ],
        "G1": [2, 0, 4, 3, 1, 2],
        "G2": [0, 0, 0, 0, 0, 0],
        "G3": [0, 0, 0, 0, 0, 0],
        "G4": [0, 26, 4, 3, 1, 1],
        "G5": [0, 0, 0, 0, 0, 0],
        "G6": [0, 0, 0, 0, 0, 0],
        "cluster index": [2, 3, 4, 5, 6, 7],
    }
    df = pd.DataFrame(data)
    path = tmp_path / "gnps_network.tsv"
    df.to_csv(path, sep="\t", index=False)
    return path


@pytest.fixture
def mock_direct_gnps(tmp_path):
    """Fixture for the direct GNPS format with '#ClusterIdx' and '#Filename'."""
    data = {
        "#ClusterIdx": [61, 61, 62, 62, 62],
        "#Filename": [
            "input_spectra/sampleA.mzXML",
            "input_spectra/file_ref1.mzXML",
            "input_spectra/sampleB.mzML",
            "input_spectra/file_ref2.mzXML",
            "input_spectra/file_ref1.mzXML",
        ],
        "#SpecIdx": [100, 200, 300, 400, 500],
        "#Scan": [111, 222, 333, 444, 555],
        "#ParentMass": [100.5, 101.1, 102.2, 103.3, 104.4],
        "#Charge": [0, 0, 1, 1, 1],
        "#RetTime": [1.1, 2.2, 3.3, 4.4, 5.5],
        "#PrecIntensity": [1000, 2000, 3000, 4000, 5000],
    }
    df = pd.DataFrame(data)
    path = tmp_path / "direct_gnps.tsv"
    df.to_csv(path, sep="\t", index=False)
    return path


@pytest.fixture
def mock_reference_metadata(tmp_path):
    """Fixture for external reference metadata with ontology columns."""
    df = pd.DataFrame(
        {
            "filename": [
                "file_ref1.mzXML",
                "file_ref2.mzXML",
                "file_ref3.mzXML",
            ],
            "ontology_terminal_leaf": ["complex", "plant", "water"],
            "sample_type_group1": ["complex", "plant", "water"],
            "sample_type_group2": ["complex", "fruit", "water"],
            "sample_type_group3": ["complex", "fruit", "water"],
            "sample_type_group4": ["complex", "fruit", "water"],
            "sample_type_group5": ["complex", "fruit", "water"],
            "sample_type_group6": ["complex", "fruit", "water"],
            "sample_name": ["complex", "fruit", "water"],
        }
    )
    path = tmp_path / "reference_metadata.csv"
    df.to_csv(path, index=False)
    return path


@pytest.fixture
def mock_sample_metadata(tmp_path):
    """Fixture for external sample metadata."""
    df = pd.DataFrame(
        {
            "filename": ["sampleA", "sampleB"],
            "group": ["G1", "G2"],  # Ensure sample names match those in GNPS
        }
    )
    path = tmp_path / "sample_metadata.csv"
    df.to_csv(path, index=False)
    return path


@pytest.fixture
def load_test_files(
    mock_network_based_gnps,
    mock_direct_gnps,
    mock_sample_metadata,
    mock_reference_metadata,
):
    def _load(fixture_type):
        if fixture_type == "gnps1":
            return (
                mock_network_based_gnps,
                mock_sample_metadata,
                mock_reference_metadata,
            )
        elif fixture_type == "gnps2":
            return (
                mock_direct_gnps,
                mock_sample_metadata,
                mock_reference_metadata,
            )
        else:
            raise ValueError(
                "Unknown fixture type. Use 'network' or 'direct'."
            )

    return _load
