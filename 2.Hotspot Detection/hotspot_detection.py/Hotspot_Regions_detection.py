# ==========================
# Packages
# ==========================
import os
import glob
import gzip
import argparse
import numpy as np
import pandas as pd 

# ==========================
# Argument Parsing
# ==========================

parser = argparse.ArgumentParser(description='Process directories and detect hotspots.')

parser.add_argument('directory_path')# Argument for the directory path containing folders for each chromosome, wich inside contains .tsv.gz file for each gene.
parser.add_argument('output_filepath')# Argument for the output file path where the results will be saved.

args = parser.parse_args()

# ==========================
# Directory Path
# ==========================
dir_path = args.directory_path

# ==========================
# Hotspot Detection Function
# ==========================
# Function to detect hotspots in a file based on codon index and fractions.
def Hotspot_detection(file_path, chromosome):
    # Extract the gene and transcript from the file name.
    file_name = os.path.basename(file_path)
    gene, transcript = file_name.split('_')[0], file_name.split('_')[1].split('.')[0]  # Expected format: gene_transcript
    
    # Read the TSV file into a DataFrame.
    with gzip.open(file_path, 'rt') as f:
        df = pd.read_csv(f, sep='\t')
    
    # Verify that the required column 'Hotspot_u_nmers002' exists.
    if 'Hotspot_u_nmers002' not in df.columns:
        return pd.DataFrame()  # Return an empty DataFrame if the column is missing.
    
    # Dictionary to store results for each gene and transcript.
    dic = {'Gene': [], 'Transcript': [], 'Gene_len': [], 'Htspt_len': [], 'Htspt_coord': [], 'Codon_GAP': [], 'Mean_fracc': [], 'Chromosome': []}
    
    # Variables to track sequences of 'YES' in the data.
    yes = []  # List to store lengths of 'YES' sequences.
    coord = []  # List to store the coordinates of 'YES' sequences.
    GAP = []  # List to store codon gaps.
    fracc = []  # List to store fraction values for calculating the mean.
    start_coord = None  # Start coordinate for 'YES' sequence.
    end_coord = None  # End coordinate for 'YES' sequence.
    yes_count = 0  # Counter for 'YES' codons.
    no_count = 0  # Counter for 'NO' codons.
    gene_len = len(df)  # Total length of the gene.
  
    
    # Loop through each row in the DataFrame.
    for idx, row in df.iterrows():
        if row['Hotspot_u_nmers002'] == 'YES':
            no_count = 0  # Reset 'NO' counter when a 'YES' is found.
            yes_count += 1  # Increment 'YES' counter.
            if start_coord is None:
                start_coord = row['CodonIndex']  # Record start coordinate for 'YES' sequence.
            fracc.append(row['Fraction_u_pnmers_002'])  # Store the fraction value.
        elif row['Hotspot_u_nmers002'] == 'NO':
            no_count += 1  # Increment 'NO' counter.
            if (yes_count != 0) or (len(yes) != 0):
                GAP.append(row['CodonIndex'])  # Record the codon index as a gap.
                fracc.append(row['Fraction_u_pnmers_002'])  # Store the fraction value.
            if (no_count == 1) and (yes_count != 0):
                yes.append(yes_count)  # Store the length of the 'YES' sequence.
                end_coord = row['CodonIndex'] - 1  # Record the end coordinate.
                coord.append(f"{start_coord}-{end_coord}")  # Save the coordinates.
                yes_count = 0  # Reset the 'YES' counter.
                start_coord = None  # Reset start coordinate.
                end_coord = None  # Reset end coordinate.
            if (no_count > 2) and (len(yes) != 0):
                # Check if any 'YES' sequence is long enough to be considered a hotspot.
                if any(x >= 5 for x in yes):
                    starts = [int(item.split('-')[0]) for item in coord if item]
                    ends = [int(item.split('-')[1]) for item in coord if item]
                    if starts and ends:
                        start = starts[0]  # Start of the hotspot.
                        end = ends[-1]  # End of the hotspot.
                        fracc_value = fracc[:-no_count]  # Exclude 'NO' fractions.
                        mean_f = np.mean(fracc_value)  # Calculate mean fraction.
                        GAPs = GAP[:-no_count]  # Exclude 'NO' gaps.
                        dic["Gene"].append(gene)
                        dic["Transcript"].append(transcript)
                        dic["Gene_len"].append(gene_len)
                        dic["Htspt_len"].append(end - start + 1)
                        dic["Htspt_coord"].append(f"{start}-{end}")
                        dic["Codon_GAP"].append(','.join(map(str, GAPs)))
                        dic["Mean_fracc"].append(mean_f)
                        dic["Chromosome"].append(chromosome)
                    # Reset tracking variables for the next sequence.
                    yes = []
                    coord = []
                    GAP = []
                    fracc = []
                    start_coord = None
                    end_coord = None
                    yes_count = 0
                    no_count = 0
                else:
                    # Reset if the sequence is not a valid hotspot.
                    yes = []
                    coord = []
                    GAP = []
                    fracc = []
                    start_coord = None
                    end_coord = None
                    yes_count = 0
                    no_count = 0

    # Final check for any remaining 'YES' sequences at the end of the file.
    if (row['CodonIndex'] == gene_len) and (yes_count != 0):
        yes.append(yes_count)
        end_coord = row['CodonIndex']
        if any(x >= 5 for x in yes):
            coord.append(f"{start_coord}-{end_coord}")
            starts = [int(item.split('-')[0]) for item in coord]
            ends = [int(item.split('-')[1]) for item in coord]
            start = starts[0]
            end = ends[-1]
            mean_f = np.mean(fracc)
            GAPs = GAP[:-no_count]
            dic["Gene"].append(gene)
            dic["Transcript"].append(transcript)
            dic["Gene_len"].append(gene_len)
            dic["Htspt_len"].append(end - start + 1)
            dic["Htspt_coord"].append(f"{start}-{end}")
            dic["Codon_GAP"].append(','.join(map(str, GAPs)))
            dic["Mean_fracc"].append(mean_f)
            dic["Chromosome"].append(chromosome)
        else:
            # Reset if the sequence is not a valid hotspot.
            yes = []
            coord = []
            GAP = []
            fracc = []
            start_coord = None
            end_coord = None
            yes_count = 0
            no_count = 0
                
    if (row['CodonIndex'] == gene_len) and (row['Hotspot_u_nmers002'] == 'NO') and (len(yes) != 0):
        if any(x >= 5 for x in yes):
            starts = [int(item.split('-')[0]) for item in coord]
            ends = [int(item.split('-')[1]) for item in coord]
            start = starts[0]
            end = ends[-1]
            mean_f = np.mean(fracc[:-no_count])
            GAPs = GAP[:-no_count]
            dic["Gene"].append(gene)
            dic["Transcript"].append(transcript)
            dic["Gene_len"].append(gene_len)
            dic["Htspt_len"].append(end - start + 1)
            dic["Htspt_coord"].append(f"{start}-{end}")
            dic["Codon_GAP"].append(','.join(map(str, GAPs)))
            dic["Mean_fracc"].append(mean_f)
            dic["Chromosome"].append(chromosome)
        else:
            # Reset if the sequence is not a valid hotspot.
            yes = []
            coord = []
            GAP = []
            fracc = []
            start_coord = None
            end_coord = None
            yes_count = 0
            no_count = 0

    # Convert the dictionary into a DataFrame and return it.
    result_df = pd.DataFrame(dic)
    return result_df

# ==========================
# Process Directory Function
# ==========================
# Function to process all files in a chromosome directory.
def process_directory(chromosome, dir_path):
    results = []
    file_paths = glob.glob(f"{dir_path}/{chromosome}/scores/*.tsv.gz")  # Get all .tsv.gz files in the directory.
    
    for file_path in file_paths:
        # Process each file and append the results.
        result = Hotspot_detection(file_path, chromosome)
        results.append(result)
    
    if results:
        # Concatenate results into a single DataFrame.
        final_result = pd.concat(results, ignore_index=True)
        return final_result
    else:
        # Return an empty DataFrame if no results.
        return pd.DataFrame()

# ==========================
# Main Script Execution
# ==========================
# Define all chromosome directories to be processed.
chromosomes = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']

all_chromosomes_results = []

# Loop through each chromosome and process the respective directory.
for chromosome in chromosomes:
    print(f"Processing directory for {chromosome}")
    chromosome_result = process_directory(chromosome, dir_path)
    if not chromosome_result.empty:
        all_chromosomes_results.append(chromosome_result)
    print(f"Finished processing directory for {chromosome}")

# Concatenate all results into a single DataFrame.
final_df = pd.concat(all_chromosomes_results, ignore_index=True)

# Sort the final DataFrame by Gene and Transcript.
final_df = final_df.sort_values(by=['Gene', 'Transcript'])

# Reorder the columns
final_df = final_df[['Chromosome', 'Gene', 'Transcript', 'Gene_len', 'Htspt_len', 'Htspt_coord', 'Codon_GAP', 'Mean_fracc']]

# Save the final DataFrame to a TSV file.
final_df.to_csv(args.output_filepath, sep='\t', index=False)
