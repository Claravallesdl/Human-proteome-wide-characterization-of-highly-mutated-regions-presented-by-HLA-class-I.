## **README (plots_distances_all_clusters.py)**

Plots and summary statistics of the distances between the HPC clusters and the distances bewteen the first/last codon of the gene and the first/last HPC cluster (tail distances).

### **Run file** 
`run_distances_all_clusters.sh`

### **Input Files** 
`'gene'_'transcript'.tsv.gz`

The input for this code consists of .tsv files, stored by chromosome in different folders. There is one file for each gene in the proteome, resulting in a total of 19,205 .tsv files.

### **Output Files**
*001/022

**`plot_distances_all_clusters_lhr*.png`**

Violin plot of the distances between the HCP clusters.

**`plot_tail_distances_lhr.png`**

Violin plot of the distances bewteen the first/last codon of the gene and the first/last HPC cluster (tail distances).

**`plot_zoom_distances_all_clusters_lhr*.png`**

Violin plot zoomed of the distances between the HPC clusters.

**`statistics__distances_between_clusters_lhr.txt`**

Statistics summary of the distances between clusters:

- Mean
- Median
- Q1
- Q3
- Min
- Max

**`statistics__tail_distances_lhr.txt`**

Statistics summary of the tail distances:

- Mean
- Median
- Q1
- Q3
- Min
- Max

### **Code Structure**

### **Functions explanation**

