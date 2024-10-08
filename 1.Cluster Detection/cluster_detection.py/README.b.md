## **README (cluster_detection.py)**

Detection of Highly Presented Codons (HPC) across all the human proteome.

### **Run file** 
`run_cluster_detection.sh`

### **Input Files** 
`'gene'_'transcript'.tsv.gz`

The input for this code consists of .tsv files. There is one file for each gene in the proteome, meaning the input will be a total of **19,205 .tsv files**.

### **Output File**

`cluster_detection_(001/002).tsv.gz

One single file with all the genes. There is a row for each cluster length for each gene. It means taht it can be more than one row per gene.

`Chromosome:` Chormosome

`Gene:` Gene name extracted from the input file name.

`Transcript:` Transcript name extracted from the input file name.

`Gene_Length:` Gene length in aa (number of rows of the file).

`Num_Clusters:` Number of clusters for an specific cluster length in the gene.

`Length_Clusters:` Number of consecutive HPC.

`Codon_coordinates:` Coordinates of the cluster (codon index of the fisrt and the last HPC that constitutes the cluster).

`Peptides:` The total number of presented peptides generated by the codons included in the cluster, calculated as the sum of the 'Peptides' column from all rows included in the cluster.

`Fraction_mean:` Mean value of the 'Fraction_u_pnmers_002' column for each cluster.

`Fraction_std:` standard deviation of the 'Fraction_u_pnmers_002' column for each cluster.

### **Code Structure**
1. **Packages**
    * os
    * glob
    * argparse
    * pandas (pd)
2. Functions definition
   * HPC_clusters
   * caclulate_metrics
   * process_directory (it apply _HPC_clusters_ and _calculate_metrics_ functions)
3. List of chr to process
4. Empty list to store results for all chr (all_chromosomes_results)
5. Process each chr directory (process_directory function)
6. Concatenate results (stored in all_chromosomes_results)
7. Save the Dataframe as TSV file

### **Functions explanation**

1. **HPC_clusters function:** This function analyzes the HPC clusters within a given DataFrame and returns a summary DataFrame containing the following columns:
    *   **Num\_Hotspots**Number of hotspot regions of each size.
    *   **Length\_Hotspots** Sizes of the hotspot regions.
    *   **Codon\_coordinates**Start and end coordinates of each hotspot region.

2. **calculate\_metrics function**: This function calculates metrics such as sum of peptides, mean fraction, and standard deviation of fraction for each cluster within the specified coordinates.
    *   **Peptides**
    *   **Fraction\_mean**
    *   **Fraction\_std**
    
3. **process\_directory function**: This function processes all .tsv files within a specified directory, calculates HPC clusters, and computes metrics for each cluster.

### **Functions flow chart**
