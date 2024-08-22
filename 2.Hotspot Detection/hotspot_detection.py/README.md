# Hotspot Detection Script

This Python script processes genomic data files organized by chromosome folders to detect hotspots regions of Highly Presented Codons (HPC). The script reads compressed TSV files, identifies significant "hotspot" regions, and outputs the results into a single TSV file.

## Requirements
- Python 3.x
- NumPy
- Pandas

## Usage

The script takes two arguments: the directory containing the chromosome subdirectories and the path where the output file will be saved.

Run the script with:

```bash
python hotspot_detection.py <directory_path> <output_filepath>
```
<directory_path>: The root directory containing subdirectories for each chromosome.
<output_filepath>: The full path where the resulting TSV file will be saved.

## Data
The input data for this script consists of compressed `.tsv.gz` files. Each file corresponds to a single gene in the proteome, totaling 19,295 files. These files are organized by chromosome into separate folders. The format of each file is a DataFrame where rows represent the codon indices of the gene, and columns contain various pieces of information. The columns relevant to this script are:

- **`Gene`**: The gene identifier.
- **`Transcript`**: The transcript identifier for the gene.
- **`Hotspot_u_nmers_002`**: Indicates whether a codon is part of a hotspot ("YES") or not ("NO").
- **`Fraction_u_pnmers_002`**: A numeric value representing the fraction.

## Functions

### `Hotspot_detection(file_path, chromosome)`
This function analyzes a single TSV file and detects hotspot regions based on specific criteria.

**Parameters:**
- **`file_path (str)`**: The path to the compressed TSV file.
- **`chromosome (str)`**: The chromosome identifier.

**Returns:**
- **`result_df (DataFrame)`**: A DataFrame containing the detected hotspots for the given file, with columns such as `Gene`, `Transcript`, `Gene_len`, `Htspt_len`, `Htspt_coord`, `Codon_GAP`, `Mean_fracc`, and `Chromosome`.

---

### `process_directory(chromosome, dir_path)`
This function processes all files in a specified chromosome directory.

**Parameters:**
- **`chromosome (str)`**: The chromosome identifier (e.g., 'chr1', 'chrX').
- **`dir_path (str)`**: The root directory containing the chromosome subdirectories.

**Returns:**
- **`final_result (DataFrame)`**: A DataFrame containing the combined hotspot detection results for all files in the chromosome's directory.

---

### Main Execution Block
The script iterates through a list of chromosomes (`chr1` to `chr22`, `chrX`, and `chrY`), processes each directory using `process_directory()`, and aggregates the results. The final combined DataFrame is then sorted and saved as a TSV file at the specified output path.

---

## How It Works

1. **Hotspot Detection**: The script reads each file, identifying sequences of "YES" in the `Hotspot_u_nmers002` column. If these sequences are at least 5 codons long, they are considered hotspots.
2. **Result Compilation**: Detected hotspots are stored in a dictionary, which is later converted into a Pandas DataFrame.
3. **Processing Directories**: The script processes all chromosome subdirectories sequentially, combining the results into a final DataFrame.
4. **Output**: The results are saved in the specified output path as a TSV file.

## Output
rows: 1 row per hotspot region detected.
The output is a TSV file with the following columns:

- **`Chromosome`**: The chromosome identifier (e.g., 'chr1', 'chrX').
- **`Gene`**: The gene name extracted from the file names.
- **`Transcript`**: The transcript name extracted from the file names.
- **`Gene_len`**: The total length of the gene.
- **`Htspt_len`**: The length of the detected hotspot region.
- **`Htspt_coord`**: The coordinates of the hotspot region.
- **`Codon_GAP`**: The gaps between codons in the hotspot region.
- **`Mean_fracc`**: The average value of the fractions for the detected hotspot.

The results are sorted by `Gene`, with potentially multiple rows per gene if more than one hotspot is detected.

