# gnps-rdd documentation

Reference data-driven metabolomics analysis package.

## Overview

This package implements the reference data-driven (RDD) metabolomics approach for analyzing untargeted tandem mass spectrometry (MS/MS) data. RDD metabolomics contextualizes experimental spectra by comparing them against curated reference datasets, allowing researchers to infer the likely origin or source of spectra even when precise chemical identities remain unknown.

The package processes molecular networking outputs from GNPS Classical Molecular Networking, matching experimental spectra against hierarchical reference databases to quantify metabolite distributions across sample groups. By leveraging hierarchical metadata (e.g., food ontologies spanning from broad categories to specific items), RDD analysis can increase MS/MS spectral usage rates from 5–6% to 15–30%.

This Python implementation provides programmatic access to RDD capabilities for advanced and large-scale analyses, with support for multivariate statistics, compositional data analysis, and flexible visualization options for exploring metabolite patterns in complex biological samples.

## Key features

- Support for GNPS1 and GNPS2 network formats
- Hierarchical ontology analysis (multi-level classification)
- Filtering by ontology level, reference type, and sample group
- Visualization tools (Matplotlib and Plotly backends)
- Statistical analysis including PCA

## Installation

From source:

```bash
git clone https://github.com/bittremieux-lab/gnps-rdd.git
cd gnps-rdd
pip install -e .
```

For development with testing tools:

```bash
pip install -e ".[dev]"
```


## Quick start

```python
from rdd import RDDCounts


# Define data source parameters
gnps_food_nist_task = "74089e95b8df41b2af7c289869dc866f"  # GNPS task ID
gnps_food_nist_path = "../data/sample_gnps_vegomn.tsv"    # Local file backup

# Attempt to load data from GNPS, fallback to local file if needed
try:
    rdd = RDDCounts(
        task_id=gnps_food_nist_task,
        gnps_2=False,                    
        sample_types="simple",           
        sample_groups=["G1"],            
        reference_groups=["G4"],         
        ontology_columns=None           
    )
    print("✓ Successfully loaded data from GNPS")
except Exception as e:
    print(f"An error occurred accessing GNPS. Will load data from file instead.")
    rdd = RDDCounts(
        gnps_network_path=gnps_food_nist_path,
        sample_types="simple",
        sample_groups=["G1"],
        reference_groups=["G4"],
        ontology_columns=None
    )

# Get counts at ontology level 3
counts = rdd.get_counts(level=3)

# Visualize
from rdd.visualization import Visualizer, MatplotlibBackend
viz = Visualizer(MatplotlibBackend())
viz.plot_reference_type_distribution(rdd, level=3)
```
## Contents

```{toctree}
:maxdepth: 2
:caption: User guide

getting_started
food_ontology
tutorials/index
```

```{toctree}
:maxdepth: 2
:caption: API reference

api/rddcounts
api/utils
api/visualization
api/analysis
```

## Indices

* {ref}`genindex`
* {ref}`modindex`
* {ref}`search`