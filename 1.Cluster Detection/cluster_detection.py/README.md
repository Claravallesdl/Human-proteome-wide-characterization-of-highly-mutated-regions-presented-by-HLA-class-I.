# HPC Cluster Detection Script

This Python script processes genomic data to analyze Highly Presented Codons (HPC) clusters in High Performance Computing (HPC) environments. It obtains the coordinates and calculates metrics for clusters of HPCs and outputs the results to a tsv file.

## Requirements

The script requires the following Python packages:

- `os`
- `glob`
- `argparse`
- `pandas`

## Usage

To run the script, use the following command:

```bash
python cluster_detection.py <directory_path> <output_filepath> <lhr>
```
## Arguments

- `<directory_path>`: Path to the directory containing chromosome folders (from 1 to 22, X and Y). Inside chromosome folders there is a .tsv file for each gene of the chromosome.
- `<output_filepath>`: Path where the output file with the final results will be saved.
- `<lhr>`: Specify lhr as `001` or `002`.

## Functions

### `HPC_clusters(dataframe)`

Identifies clusters of Highly Presented Codons (HPC) and their sizes from the input DataFrame.

**Input:**

- `dataframe (pd.DataFrame)`: DataFrame for each gene.

**Returns:**

- `pd.DataFrame`: DataFrame summarizing the number and coordinates of HPC clusters.

**Explanation:**

This function iterates through the rows of the input DataFrame to identify clusters of Highly Presented Codons (HPC). It determines whether each codon is considered highly presented based on the presence of 'YES' or 'NO' in the `Hotspot_u_nmers_lhr` column. A 'YES' indicates a highly presented codon, while a 'NO' indicates it is not highly presented. It tracks the size of each cluster and its start and end positions. The function then creates a DataFrame summarizing the number of clusters, their sizes, and their coordinates.

### `calculate_metrics(cluster_coords, data)`

Calculates various metrics for clusters based on provided coordinates.

**Input:**

- `cluster_coords (str)`: Comma-separated string of cluster coordinates.
- `data (pd.DataFrame)`: DataFrame containing genomic data.

**Returns:**

- `dict`: Dictionary containing calculated metrics as strings joined by commas.

**Explanation:**

This function calculates metrics for each cluster based on its coordinates. It computes the sum of peptides, the mean, and the standard deviation of the `Fraction_u_pnmers` values within each cluster. The results are returned as a dictionary where values are strings joined by commas.

### `process_directory(chromosome, dir_path)`

Processes all files in the specified directory for a given chromosome and computes cluster metrics.

**Input:**

- `chromosome (str)`: Chromosome identifier (e.g., 'chr1', 'chrX').
- `dir_path (str)`: Path to the directory containing chromosome folders.

**Returns:**

- `pd.DataFrame`: DataFrame containing processed results for the chromosome.

**Explanation:**

This function processes all TSV files within a specified directory for a given chromosome. It reads each file into a DataFrame, computes HPC clusters, calculates various metrics for these clusters, and then combines the results.

### Main Execution

This section of the code was designed to automate the processing of data for all chromosome folders. Instead of manually processing each chromosome separately, this loop iterates through a list of chromosomes and processes each one using the `process_directory` function.

#### Steps:

1. **List of Chromosomes**: The script creates a list of chromosome identifiers, including `chr1` through `chr22`, as well as `chrX` and `chrY`.

    ```python
    chromosomes = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']
    ```

2. **Initialize Results Storage**: An empty list `all_chromosomes_results` is initialized to store the results for each chromosome.

    ```python
    all_chromosomes_results = []
    ```

3. **Process Each Chromosome**: The script loops through each chromosome in the list, calling the `process_directory` function for each one. The processed results for each chromosome are appended to the `all_chromosomes_results` list.

    ```python
    for chromosome in chromosomes:
        chromosome_result = process_directory(chromosome, dir_path)
        all_chromosomes_results.append(chromosome_result)
    ```

4. **Concatenate Results**: After processing all chromosomes, the script concatenates the results for each chromosome into a single DataFrame. This allows the data for all chromosomes to be stored together.

    ```python
    final_result = pd.concat(all_chromosomes_results, ignore_index=True)
    ```

## Output Description

The final TSV file produced by the script contains the following columns:

- **Chromosome**: The chromosome identifier (e.g., 'chr1', 'chrX').
- **Gene**: The gene name extracted from the file names.
- **Transcript**: The transcript name extracted from the file names.
- **Gene_Length**: Length of the gene (number of rows in the input data).
- **Num_Clusters**: Number of clusters of each size.
- **Length_Clusters**: Sizes of the clusters.
- **Codon_coordinates**: Coordinates of each cluster.
- **Peptides**: Total number of presented peptides for each cluster.
- **Fraction_mean**: Mean of `Fraction_u_pnmers_001/002` for each cluster.
- **Fraction_std**: Standard deviation of `Fraction_u_pnmers_001/002` for each cluster.

The results are sorted by **Gene** and then by **Length_Clusters** within each Gene.
