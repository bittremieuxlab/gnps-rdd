# GNPS Reference Data-Driven (RDD) Metabolomics Package

A Python package for reference data-driven analysis of untargeted MS/MS metabolomics data using GNPS Classical Molecular Networking outputs.

[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

## Overview

RDD metabolomics contextualizes experimental MS/MS spectra by comparing them against curated reference datasets organized in hierarchical ontologies. Rather than identifying exact molecular structures, this approach infers the likely origin or source of spectra by matching them to reference categories through GNPS molecular networking, even when precise chemical identities remain unknown.

This method increases MS/MS spectral usage rates from 5-6% to 15-30%—a 3- to 5-fold improvement in data interpretability—enabling researchers to extract biological meaning from previously unannotated spectra.

## Key Features

Generate RDD counts from GNPS Classical Molecular Networking outputs with hierarchical reference matching across customizable ontology levels. The package provides flexible filtering and aggregation by ontology level, reference category, or sample group, along with a rich visualization suite supporting both Matplotlib and Plotly backends. Statistical analysis tools include PCA with compositional data support, enabling creation of publication-ready and interactive figures.

## Installation
```bash
pip install gnps-rdd
```

Or from source:
```bash
git clone https://github.com/bittremieux-lab/gnps-rdd.git
cd gnps-rdd
pip install -e .
```

## Workflow

1. **Prepare GNPS network**: Upload experimental samples and reference data with hierarchical metadata to GNPS, then run Classical Molecular Networking
2. **Load and calculate**: Use the package to load GNPS outputs and calculate RDD counts across ontology levels
3. **Analyze and visualize**: Filter, aggregate, and visualize results using built-in statistical and plotting tools

See the [GNPS RDD tutorial](https://ccms-ucsd.github.io/GNPSDocumentation/tutorials/rdd/) for detailed networking setup instructions.

## Use Cases

**Food metabolomics**: Reconstruct dietary patterns from plasma or stool samples by matching against comprehensive food reference databases (e.g., Global FoodOmics Project with 3,500+ food items across 6 hierarchical levels).

**Environmental metabolomics**: Track contamination sources and environmental exposures using location-based or source-type hierarchies.

**Microbial metabolomics**: Link metabolites to taxonomic origins through phylogenetic hierarchies (domain → phylum → class → genus → species).

## Documentation

- **Full documentation**: [https://gnps-rdd.readthedocs.io](https://gnps-rdd.readthedocs.io)
- **GNPS tutorial**: [https://ccms-ucsd.github.io/GNPSDocumentation/tutorials/rdd/](https://ccms-ucsd.github.io/GNPSDocumentation/tutorials/rdd/)
- **Web application**: [https://gnps-rdd.streamlit.app/](https://gnps-rdd.streamlit.app/)
- **Example notebooks**: [rdd/notebooks](rdd/notebooks)

## Citation

If you use this package in your research, please cite:
```
Alejandro Mendoza Cantu, Julia M. Gauglitz, Wout Bittremieux. "An open-source platform for reference data-driven analysis of untargeted metabolomics" In preparation (2025).
```

## License

Apache License 2.0 - see [LICENSE](LICENSE) file for details.

## Support

- **Issues**: [GitHub Issues](https://github.com/bittremieux-lab/gnps-rdd/issues)
