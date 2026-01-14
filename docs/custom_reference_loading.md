# Custom reference data and associated metadata

Reference data-driven (RDD) metabolomics is a flexible method that can be used in domains other than foodomics. The main factor for running an effective RDD analysis is the creation of a reference dataset with its associated curated and rich metadata. We will briefly describe the main components needed to build this reference dataset and how to load it when running an analysis with the gnps-rdd Python package. This tutorial does not delve into the specifics of experimental data acquisition, as this depends on the specific research question and the experimental setup of the researcher. It is important to note that high-resolution MS/MS data is needed for the construction of the references.

## Reference data

The reference data must be high-resolution MS/MS data that properly covers a reference space linked to a specific research question. An example is the Global Foodomics Project dataset, used to validate the use of RDD in the field of foodomics and as the default reference for the gnps-rdd Python package. The 3,527 food references included in it allowed the RDD methodology to gain insights into questions such as: Can individuals who follow an omnivore or vegan diet be separated using the information from the spectral matches with this reference dataset? Can RDD detect specific food categories used in the diet for a specific study?

## Minimum metadata requirements

The minimum metadata requirements for the reference dataset are just the filename associated with each MS/MS run and a sample name that identifies it. This will allow RDD to build a counts table at the most specific level, although, logically, there is no possibility of doing ontology aggregation analysis.

| filename                    | sample_name  |
| --------------------------- | ------------ |
| 72500_1x_BG3_01_17471.mzXML | 11442.G72500 |
| 72504_1x_BH8_01_17502.mzXML | 11442.G72504 |
| 72507_1x_BH9_01_17503.mzXML | 11442.G72507 |

## Recommended metadata requirements

The addition of an ontology classification to the reference dataset enhances the analysis that can be done with RDD. The addition of different ontology levels enables visualization and pattern identification across different categories and their aggregations.

| filename                    | sample_name  | sample_type_group1 | sample_type_group2 | sample_type_group3 | sample_type_group4 | sample_type_group5 | sample_type_group6 |
| --------------------------- | ------------ | ------------------ | ------------------ | ------------------ | ------------------ | ------------------ | ------------------ |
| 72500_1x_BG3_01_17471.mzXML | 11442.G72500 | plant              | fruit              | fleshy fruit       | pome               | apple              | apple              |
| 72504_1x_BH8_01_17502.mzXML | 11442.G72504 | plant              | fruit              | fleshy fruit       | multifruit         | fig                | turkish fig        |
| 72507_1x_BH9_01_17503.mzXML | 11442.G72507 | algae              | algae              | seaweed            | seaweed            | dulse              | dulse              |

## Metadata format

The metadata file can be uploaded from the user's local directory and should follow the structure of the examples shown above. Users can select custom column names, which can be specified when calling the RDD class. The argument `ontology_columns` accepts an array with the column names corresponding to the ontology classification. However, the columns `filename` and `sample_name` must exist in the metadata. The file can have the .txt, .csv, or .tsv extensions.

The `levels` parameter is optional. When not specified, the number of ontology levels is automatically determined from:
1. The length of `ontology_columns` if custom columns are provided, or
2. The number of `sample_type_groupX` columns in the default metadata

## Example: Loading custom metadata
```python
# Levels will be automatically set to 3 based on ontology_columns length
rdd_counts = RDDCounts(
    gnps_network_path=network_path,
    external_reference_metadata=metadata_path,
    sample_groups=["G1"],
    reference_groups=["G4"],
    ontology_columns=["level_1", "level_2", "level_3"]
)

# Or explicitly specify levels (must not exceed available ontology columns)
rdd_counts = RDDCounts(
    gnps_network_path=network_path,
    external_reference_metadata=metadata_path,
    sample_groups=["G1"],
    reference_groups=["G4"],
    ontology_columns=["level_1", "level_2", "level_3"],
    levels=2  # Only use first 2 levels
)
```
## Example: Accesing custom metadata

```python
rdd_counts.reference_metadata
```
