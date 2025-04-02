import pytest
import pandas as pd
from io import StringIO
import sys
import os
# Add the path to the 'rdd' directory to the system path
project_path = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "rdd")
)
sys.path.append(project_path)
from rdd.RDDcounts import RDDCounts
from rdd.utils import _validate_groups, _load_RDD_metadata, calculate_proportions
import pkg_resources


def test_validate_groups_invalid():
    gnps = pd.DataFrame({"DefaultGroups": ["G1", "G2"]})
    with pytest.raises(ValueError, match="invalid"):
        _validate_groups(gnps, groups_included=["G3"])


def test_calculate_proportions_raises_on_multiple_levels():
    df = pd.DataFrame({
        "filename": ["f1", "f2"],
        "reference_type": ["a", "b"],
        "count": [1, 2],
        "level": [1, 2],
        "group": ["G1", "G1"]
    })
    with pytest.raises(ValueError, match="Multiple levels found"):
        calculate_proportions(df)


def test_load_rdd_metadata_pkg_resources(monkeypatch):
    def fake_open_text_fail(*args, **kwargs):
        raise ModuleNotFoundError

    def fake_resource_stream(*args, **kwargs):
        return StringIO("filename\tsample_type\nfile1\tcomplex")

    monkeypatch.setattr("rdd.utils.resources.open_text", fake_open_text_fail)
    monkeypatch.setattr(pkg_resources, "resource_stream", fake_resource_stream)

    df = _load_RDD_metadata()
    assert not df.empty
    assert "filename" in df.columns



