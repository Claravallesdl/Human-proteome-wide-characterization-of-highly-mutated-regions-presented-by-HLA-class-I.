## Overview

This script is designed to analyze and compare disorder regions and hotspot regions within genes. It processes input data to extract the coordinates of each region, the codon indices of these regions and then compares these codon indices to identify overlaps.

The script performs the following steps:
1. **Extract Coordinates**: It extracts the coordinates of disorder and hotspot regions from the input data using regular expressions.
2. **Create Dictionaries**: It creates dictionaries where each gene is associated with sets of codon indices. These indices represent the codons that fall within the specified disorder and hotspot regions.
3. **Compare Codon Indices**: It compares the sets of codon indices for each gene to determine overlaps. Specifically, it identifies codons that are both within disordered regions and hotspot regions.

## Requirements

- Python 3.x
- `numpy`
- `pandas`
- `statsmodels`

## Data

The script expects two DataFrames:

1. **`disorder`**: Contains the coordinates of disorder regions for each gene.
The DataFrame results from merging the UniProt ID - gene annotation dataset with the `UP000005640_9606_region.bed.gz` file, available for download from [UniProt's FTP server](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/genome_annotation_tracks/UP000005640_9606_beds/), and an additional `regions.bed` file.
Relevant columns:
- `Gene`: Gene name 
- `Annot_description`: coordinate information as strings.

2. **`df_HR`**: Contains the coordinates of hotspot regions for each gene.
Relevant columns:
- `Gene`: Gene name 
- `Htspt_coord`: coordinate information as strings.

## Functions

### `extract_coordinates_disorder(coord_str)`

Extracts a list of coordinates from a disorder region string.

- **Input**: A string with format like `"P38-R80; Disordered"`.
- **Output**: A list of integers representing the coordinates from `start` to `end` (inclusive).

### `extract_coordinates_hr(coord_str)`

Extracts a list of coordinates from a hotspot region string.

- **Input**: A string with format like `"38-80"`.
- **Output**: A list of integers representing the coordinates from `start` to `end` (inclusive).

### `create_coordinate_dict(df, coord_col, coord_extractor)`

Creates a dictionary of coordinates for each gene based on a specified DataFrame and coordinate extraction function.

- **Input**:
  - `df`: DataFrame containing gene data.
  - `coord_col`: Name of the column containing coordinate strings.
  - `coord_extractor`: Function to extract coordinates from a string.
- **Output**: A dictionary with gene identifiers as keys and sets of codon indices as values.

### `compare_codons(disorder_dict, hr_dict)`

Compares disordered codons with hotspot region codons and summarizes the findings.

- **Input**:
  - `disorder_dict`: Dictionary of disordered codons by gene.
  - `hr_dict`: Dictionary of hotspot region codons by gene.
- **Output**: A DataFrame with columns:
  - `Gene`: Gene identifier.
  - `Disordered_codons`: Number of disordered codons.
  - `HR_codons`: Number of hotspot region codons.
  - `Dis_HR_codons`: Number of overlapping codons.

## Statistical Analysis

This section performs statistical analysis to evaluate the distribution of disordered codons relative to hotspot regions (HR).

1. **Calculations**:
   - **Disordered Codons Outside HR**: Computes the number of disordered codons not overlapping with hotspot regions.
   - **Disordered Codons Inside HR**: Computes the number of disordered codons that overlap with hotspot regions.
   - **Total Codons Inside HR**: Computes the total number of codons within hotspot regions.
   - **Total Codons Outside HR**: Computes the total number of codons outside hotspot regions.

2. **Proportions and Fractions**:
   - Calculates the proportion and fraction of disordered codons inside and outside hotspot regions.
   - Converts these proportions into percentages for better interpretability.

3. **Statistical Testing**:
   - Performs a Z-test to determine if there is a statistically significant difference in the proportion of disordered codons inside versus outside hotspot regions.
   - Calculates a 95% confidence interval for the difference in proportions.

The results include the Z-statistic, p-value, observed difference in proportions, and the confidence interval for the difference, providing insight into the statistical significance of the observed distribution of disordered codons.
