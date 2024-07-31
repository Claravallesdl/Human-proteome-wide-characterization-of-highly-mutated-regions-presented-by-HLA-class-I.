
import os
import glob
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='New file with ')

parser.add_argument('directory_path')   # Directory of all the files for each chromosome folder
parser.add_argument('output_filepath')  # Output path
parser.add_argument('lhr', help='lhr 001 o 002')

args = parser.parse_args()

dir_path = args.directory_path

def HPC_clusters(dataframe):
    cluster_sizes = {}  # Dictionary to store the sizes of clusters
    current_cluster_size = 0  # Initialize the size counter for the current cluster
    start = None  # Variable to track the starting codon index of a cluster
    
    for index, row in dataframe.iterrows():  # Iterate through each row of the input DataFrame
        if row['Hotspot_u_nmers'+str(args.lhr)] == 'YES':  # Check if the current row is a Highly Presented Codon (HPC)
            if current_cluster_size == 0:   # Start of a new cluster
                start = row['CodonIndex']  # Set the start index of the cluster
            current_cluster_size += 1  # Increment the size of the current cluster
        else:  # 'NO' encountered (Not a HPC)
            if current_cluster_size > 0:  # Check if there was an ongoing cluster
                end = row['CodonIndex'] - 1  # Determine the end index of the cluster
                length = current_cluster_size  # Length of the cluster
                
                if current_cluster_size in cluster_sizes:
                    # Append the start-end coordinate to the list of clusters of this size
                    cluster_sizes[current_cluster_size].append(f"{start}-{end}")
                else:
                    # Create a new list of coordinates for this cluster size
                    cluster_sizes[current_cluster_size] = [f"{start}-{end}"]
                current_cluster_size = 0  # Reset the cluster size counter
    
    # Final check for the last cluster if the last row of the Dataframe ends with 'YES'
    if current_cluster_size > 0:
        end = row['CodonIndex']  # Finalize the end index of the last cluster
        length = end - start + 1  # Calculate the length of the last cluster
        if current_cluster_size in cluster_sizes:
            # Append the start-end coordinate to the list of clusters of this size
            cluster_sizes[current_cluster_size].append(f"{start}-{end}")
        else:
            # Create a new list of coordinates for this cluster size
            cluster_sizes[current_cluster_size] = [f"{start}-{end}"]
    
    # Create a DataFrame from the dictionary of cluster sizes
    result_df = pd.DataFrame({
        'Num_Clusters': [len(v) for v in cluster_sizes.values()],           # Count of clusters for each size
        'Length_Clusters': list(cluster_sizes.keys()),                      # Sizes of the clusters
        'Codon_coordinates': [",".join(v) for v in cluster_sizes.values()]  # Coordinates of each cluster
    })
    
    return result_df  # Return the DataFrame summarizing clusters


def calculate_metrics(cluster_coords, data):
    # Initialize lists to store calculated metrics
    peptides_sum = []      # List to store the sum of peptides for each cluster
    fraction_means = []    # List to store the mean of Fraction_u_pnmers_001/002 for each cluster
    fraction_stds = []     # List to store the standard deviation of Fraction_u_pnmers_001/002 for each cluster
    
    # Iterate over each cluster coordinate provided
    for coord in cluster_coords.split(','):
        # Split the coordinate string into start and end positions
        start, end = map(int, coord.replace(',', '').split('-'))  # Replace commas and then split
        
        # Filter data within the specified cluster
        cluster_data = data[(data['CodonIndex'] >= start) & (data['CodonIndex'] <= end)]
        
        # Calculate sum of peptides for the cluster and append to the list
        cluster_peptides_sum = cluster_data['Peptides'].sum()
        peptides_sum.append(cluster_peptides_sum)
        
        # Calculate mean of Fraction_u_pnmers_001/002 for the cluster and append to the list
        fraction_means.append(cluster_data['Fraction_u_pnmers'+ str(args.lhr)].mean())
        
        # Calculate standard deviation of Fraction_u_pnmers_001/002 for the cluster and append to the list
        fraction_stds.append(cluster_data['Fraction_u_pnmers'+ str(args.lhr)].std())
    
    # Return a dictionary containing calculated metrics as strings joined by commas
    return {
        'Peptides': ','.join(map(str, peptides_sum)),               # Join peptide sums as a string
        'Fraction_mean': ','.join(map(str, fraction_means)),        # Join mean values as a string
        'Fraction_std': ','.join(map(str, fraction_stds))           # Join standard deviations as a string
    }


def process_directory(chromosome, dir_path):
    results = []
    
    # Obtain the list of all .tsv files in the specified directory
    file_paths = glob.glob(f"{dir_path}/{chromosome}/scores/*.tsv.gz")
    
    for file_path in file_paths:
        # Extract Gene and Transcript information from the file name
        file_name = file_path.split('/')[-1]
        file_name_parts = file_name.split('_')
        Gene = file_name_parts[0].split('\\')[-1]  # Extract Gene from the file name (first part before '_')
        Transcript = file_name_parts[1].split('.')[0]  # Extract Transcript from the file name (first part before '.')
        
        # Read input file into a DataFrame
        data = pd.read_csv(file_path, sep='\t')
        
        # Calculate cluster sizes and counts for the input DataFrame
        result = HPC_clusters(data)
        
        # Add Gene and Transcript columns to the result DataFrame
        result['Gene'] = Gene
        result['Transcript'] = Transcript
        result['Chromosome'] = chromosome  # Add Chromosome column
        
        # Calculate metrics (Peptides, Fraction_mean, Fraction_std) for each cluster
        result_metrics = result['Codon_coordinates'].apply(lambda x: calculate_metrics(x, data))
        metrics_df = pd.DataFrame(result_metrics.tolist(), index=result.index)
        
        # Concatenate metrics DataFrame with the result DataFrame
        final_result = pd.concat([result, metrics_df], axis=1)
        
        # Calculate gene length
        gene_length = len(data)
        
        # Add Gene_Length column
        final_result['Gene_Length'] = gene_length
        
        # Append the processed result to the list of results
        results.append(final_result)
    
    # Concatenate all results into a single DataFrame
    final_result_concat = pd.concat(results, ignore_index=True)
    
    # Sort the final DataFrame by 'Gene' and then by 'Length_Clusters' within each 'Gene'
    final_result_concat.sort_values(by=['Gene', 'Length_Clusters'], ascending=[True, True], inplace=True)
    
    # Reorder the columns as specified
    final_result_concat = final_result_concat[['Chromosome', 'Gene', 'Transcript', 'Gene_Length', 'Num_Clusters', 'Length_Clusters', 'Codon_coordinates', 
                                               'Peptides', 'Fraction_mean', 'Fraction_std']]
    
    return final_result_concat


# List of chromosomes to process
chromosomes = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']

# Initialize an empty list to store the results for all chromosomes
all_chromosomes_results = []

# Process each chromosome directory
for chromosome in chromosomes:
    chromosome_result = process_directory(chromosome, dir_path)
    all_chromosomes_results.append(chromosome_result)

# Concatenate results for all chromosomes into a single DataFrame
final_result = pd.concat(all_chromosomes_results, ignore_index=True)

# Save the final DataFrame to a TSV file
final_result.to_csv(args.output_filepath, sep='\t', index=False)
