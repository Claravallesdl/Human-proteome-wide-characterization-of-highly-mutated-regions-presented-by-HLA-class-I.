# README (hotspots_detection)

This script detect hotspot regions within all the genes that the repsective folder contains.

  

### **Description**

**run** run\_chr1\_001/002

  

**input** Gene files (Carpeta Ana Chr/Socres/). The input of this script are the diferent .tsv files that an specific folder of chise contains.

  

**output** chr1\_hotspots001.tsv.gz. The output of this script is a unique file for all the genes analyzed.

Each row represents a unique combination of gene, transcript, and hotspot region length .

  

The columns, or variables, included in the output are as follows:

  

**Gene** Gene name extracted from the input file name.

**Transcript** Transcript name extracted from the input file name.

**Gene\_Length**: Indicates the length of the gene, which is calculated by determining the number of rows in the input data corresponding to that gene.

**Num\_Hotspots**Represents the total number of hotspot regions found within the gene's transcript. This value is obtained by aggregating the counts of hotspot regions for each length category.

**Length\_Hotspots** Refers to the size (length) of each hotspot region found within the gene's transcript.

**Codon\_coordinates** Provides the start and end codon indices of each hotspot region. These coordinates are obtained by iterating through the data and identifying consecutive rows marked as hotspots.

**Peptides** Represents the sum of peptides for each hotspot region. It is calculated by summing the 'Peptides' column values within each hotspot region.

**Fraction\_mean** Represents the mean value of the 'Fraction\_u\_pnmers\_002' column for each hotspot region. This metric indicates the average fraction of u-polymerase sites for the region.

**Fraction\_std** Denotes the standard deviation of the 'Fraction\_u\_pnmers\_002' column for each hotspot region. This metric provides information about the variability or dispersion of u-polymerase sites within the region.

### **Functions**

1. **Hotspot\_regions function**: This function analyzes hotspot regions within a given DataFrame and returns a summary DataFrame containing the following columns:
    *   **Num\_Hotspots**Number of hotspot regions of each size.
    *   **Length\_Hotspots** Sizes of the hotspot regions.
    *   **Codon\_coordinates**Start and end coordinates of each hotspot region.
2. **calculate\_metrics function**: This function calculates metrics such as sum of peptides, mean fraction, and standard deviation of fraction for each region within the specified coordinates.
    *   **Peptides**
    *   **Fraction\_mean**
    *   **Fraction\_std**
3. **process\_directory function**: This function processes all .tsv files within a specified directory, calculates hotspot regions, and computes metrics for each region. The processed data is stored in a DataFrame and saved to a TSV file.