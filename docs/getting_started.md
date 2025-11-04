# Getting started

## What is reference data-driven metabolomics?

Reference data-driven (RDD) metabolomics contextualizes experimental MS/MS spectra by comparing them against curated reference datasets organized in hierarchical ontologies. Rather than identifying exact molecular structures, this approach infers the likely origin or source of spectra by matching them to reference categories at different classification levels, even when precise chemical identities remain unknown.

## How it works

### Input data

Three components are required:

- **GNPS network data**: Molecular networking results from GNPS Classical Molecular Networking (task ID or downloaded files)
- **Sample group assignment**: Experimental data assigned to sample groups during GNPS networking
- **Reference metadata**: Hierarchical classification linking reference spectra to their known origins (structure depends on your reference dataset)

### Matching process
```
GNPS Molecular Network
├─ Experimental spectra (sample group)
└─ Reference spectra (reference group)
    │
    └─ Clustered together based on spectral similarity
        │
        └─ Experimental spectra inherit reference metadata via cluster membership
            │
            └─ Hierarchical metadata propagation:
                ├─ Level 0: Specific files
                ├─ Level 1: Major class
                ├─ Level 2: Intermediate class
                ├─ Level 3: Subclass
                └─ Levels 4+: More specific levels (if available)
```

When experimental and reference spectra cluster together based on spectral similarity, the experimental spectrum inherits the hierarchical metadata from the reference spectrum.

### Count calculation

For each experimental MS run and ontology level, the method calculates:

- **RDD counts**: Number of spectra assigned to each reference category through cluster membership
- **Proportions**: Relative abundance of each category normalized by total assigned spectra

This creates an RDD count table in long format linking each experimental run to specific metadata categories across multiple hierarchical levels.

### Analysis workflow

1. Generate molecular network on GNPS with experimental data in sample groups and reference data in reference groups
2. Load GNPS networking outputs and hierarchical reference metadata
3. Calculate RDD counts at all ontology levels by matching cluster assignments to metadata
4. Filter and aggregate counts by level, reference category, or sample group
5. Visualize distributions and hierarchical relationships
6. Perform multivariate statistical analysis (PCA with optional centered log-ratio transformation)

## Basic example
```python
from rdd import RDDCounts
from rdd.visualization import Visualizer, MatplotlibBackend

# Load GNPS network data and calculate RDD counts
rdd = RDDCounts(
    task_id="your_gnps_task_id",  # GNPS networking job ID
    sample_types="reference_metadata.tsv",  # Hierarchical metadata file
    sample_groups=["G1"],  # Experimental sample group(s)
    reference_groups=["G4"]  # Reference group(s)
)

# Get RDD counts at level 2 (intermediate classification level)
counts = rdd.get_counts(level=2)
print(counts.head())

# Filter for top 10 most abundant reference categories
filtered = rdd.filter_counts(level=2, top_n=10, top_n_method="total")

# Visualize distribution across samples
viz = Visualizer(MatplotlibBackend())
viz.plot_reference_type_distribution(rdd, level=2, top_n=10)

# Compare proportions between sample groups
viz.box_plot_RDD_proportions(rdd, level=2, group_by=True)
```

## Key concepts

### Hierarchical metadata levels

The reference metadata contains multiple hierarchical levels. The number and meaning of levels depends on your reference dataset:

**Example: Food reference data (Global FoodOmics Project)**
- **Level 0**: Broadest category (e.g., "Plant-based foods", "Animal-based foods")
- **Levels 1-3**: Intermediate classifications (e.g., Fruits → Citrus → Oranges)
- **Levels 4-6**: Most specific classifications (e.g., specific cultivars or preparations)


**Example: Microbial metabolite reference**
- **Level 0**: Domain (e.g., "Bacteria", "Fungi")
- **Levels 1-3**: Taxonomic ranks (Phylum → Class → Order)
- **Level 4+**: Species and strain levels



### Sample vs reference groups

- **Sample groups**: Experimental data being analyzed (assigned during GNPS networking)
- **Reference groups**: Curated reference datasets with known origins (assigned during GNPS networking)

Experimental spectra that cluster with reference spectra inherit the reference metadata through their cluster membership.

### Counts vs proportions

- **RDD counts**: Absolute number of spectra assigned to each reference category
- **Proportions**: Percentage of total assigned spectra (recommended for comparing samples with different spectral counts)

For compositional data analysis, apply centered log-ratio transformation before PCA.

## Next steps

- See the [basic usage tutorial](tutorials/example_notebook.md) (Jupyter notebook walkthrough)
- Review the API documentation:
  - [RDDCounts](api/rddcounts.md) - Main analysis class
  - [Utilities](api/utils.md) - Data processing functions
  - [Visualization](api/visualization.md) - Plotting functions
  - [Analysis](api/analysis.md) - Statistical methods
- Understand [statistical analysis](tutorials/statistical_analysis.md) options
- See the [food ontology example](examples/food_ontology.md) (Global FoodOmics use case)