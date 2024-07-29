**Objective**

Detection of Highly Presented Codons across all the human proteome.

**Run file** `run_cluster_detection.sh`

**Input Files** `gene_transcript.tsv.gz`

The input for this code consists of .tsv files. There is one file for each gene in the proteome, meaning the input will be a total of '19,205 .tsv files'.

**Output File** `cluster_detection.tsv.gz`

One single file with all the genes. There is a row for each cluster length for each gene. It means taht it can be more than one row per gene.

`Chromosome: ` Chormosome

`Gene: ` Gene name

`Transcript: ` Transcript name

`Gene_Length: ` Gene length in aa

`Num_Clusters: ` Number of clusters for an specific length in the gene

`Length_Clusters: ` Length of the cluster in aa


`Codon_coordinates: ` Coordinates of the clusters

`Peptides: ` pepitdes 

`Fraction_mean`
`Fraction_std`
