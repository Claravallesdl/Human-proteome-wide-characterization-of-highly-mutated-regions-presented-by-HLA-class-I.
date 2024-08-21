# HPC Cluster Analysis Script

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
python script_name.py <directory_path> <output_filepath> <lhr>

## Arguments

- `<directory_path>`: Path to the directory containing chromosome folders with TSV files.
- `<output_filepath>`: Path where the output file with the final results will be saved.
- `<lhr>`: Specify lhr as `001` or `002`.
