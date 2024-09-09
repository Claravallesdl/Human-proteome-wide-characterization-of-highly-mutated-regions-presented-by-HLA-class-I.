#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# ==========================
# Packages
# ==========================
import re # modue for regular expressions
import numpy as np
import pandas as pd
from statsmodels.stats.proportion import proportions_ztest, proportion_confint


# ==========================
# Data 
# ==========================

# disorder: df containing the coordinates of disorder regions of each gene
# df_HR: df conatining the coordinates of hotspot regions of each gene


# =================================
# Coordinate Extraction Functions
# =================================

# Function to extract coordinates from a disorder region string
def extract_coordinates_disorder(coord_str):
    if isinstance(coord_str, str):
        # Define the regular expression pattern to match the coordinate range
        pattern = r'([A-Za-z])(\d+)-([A-Za-z])(\d+)'
        
        # Search for the pattern in the coordinate string
        match = re.search(pattern, coord_str)
        
        if match:
            # Extract start and end coordinates
            start_letter, start_number, end_letter, end_number = match.groups()
            
            # Convert to integers
            start_coord = int(start_number)
            end_coord = int(end_number)
            
            # Return a list of coordinates
            return list(range(start_coord, end_coord + 1))
    # If not a string or no match, return an empty list
    return []

# Function to extract coordinates from a hotspot region string
def extract_coordinates_hr(coord_str):
    # Check if the input is a string
    if isinstance(coord_str, str):
        # Split the string by '-' and convert the start and end positions to integers
        start, end = map(int, coord_str.split('-'))
        # Return a list of positions from start to end (inclusive)
        return list(range(start, end + 1))
    # Return an empty list if the input is not a string
    return []

# Function to create a dictionary of coordinates for each gene
def create_coordinate_dict(df, coord_col, coord_extractor):
    coord_dict = {}  # Initialize an empty dictionary to store gene coordinates
    
    # Iterate over each row in the DataFrame
    for _, row in df.iterrows():
        gene = row['Gene']  # Extract the gene identifier from the current row
        coordinates = coord_extractor(row[coord_col])  # Extract the coordinates using the provided extractor function
        
        # Check if the gene is already in the dictionary
        if gene not in coord_dict:
            # If not, add the gene to the dictionary with the extracted coordinates
            coord_dict[gene] = set(coordinates)
        else:
            # If the gene is already present, update the set of coordinates for this gene
            coord_dict[gene].update(coordinates)
    
    return coord_dict  # Return the dictionary containing genes and their corresponding coordinates

# Apply the create_coordinate_dict function to the 'disorder' DataFrame using 'Annot_description' column
# and the extract_coordinates_disorder function to extract coordinates
disorder_dict = create_coordinate_dict(disorder, 'Annot_description', extract_coordinates_disorder)

# Apply create_coordinate_dict function to the 'df_HR' DataFrame using 'Htspt_coord' column
# and the extract_coordinates_hr function to extract coordinates
hr_dict = create_coordinate_dict(df_HR, 'Htspt_coord', extract_coordinates_hr)


# ==========================
# Comparison of Codons
# ==========================

# Function to create a DataFrame comparing disordered codons and hotspot region codons
def compare_codons(disorder_dict, hr_dict):
    results = []  # Initialize an empty list to store the results for each gene
    
    # Get the union of all gene identifiers from both disorder and hotspot dictionaries
    all_genes = set(disorder_dict.keys()).union(set(hr_dict.keys()))
    
    # Iterate over each gene in the union of both dictionaries
    for gene in all_genes:
        # Get the set of disordered codons for the gene, or an empty set if the gene is not in the disorder dictionary
        disordered_codons = disorder_dict.get(gene, set())
        # Get the set of hotspot region codons for the gene, or an empty set if the gene is not in the hotspot dictionary
        hr_codons = hr_dict.get(gene, set())
        # Find the intersection of disordered codons and hotspot region codons for the gene
        overlapping_codons = disordered_codons.intersection(hr_codons)
        
        # Append the results for the gene to the results list as a dictionary
        results.append({
            'Gene': gene,                      # Gene identifier
            'Disordered_codons': len(disordered_codons),  # Number of disordered codons
            'HR_codons': len(hr_codons),      # Number of hotspot region codons
            'Dis_HR_codons': len(overlapping_codons)  # Number of disordered codons within hotspot regions
        })
    
    # Convert the list of results dictionaries into a pandas DataFrame and return it
    return pd.DataFrame(results)

# Apply the compare_codons function to the disorder_dict and hr_dict dictionaries
disordered_codons = compare_codons(disorder_dict, hr_dict)


# ==========================
# Statistical Analysis
# ==========================

# Calculate the number of disordered codons outside hotspot regions (HR)
disordered_out_HR = disordered_codons['Disordered_codons'].sum() - disordered_codons['Dis_HR_codons'].sum()
# Calculate the number of disordered codons inside hotspot regions (HR)
disordered_in_HR = disordered_codons['Dis_HR_codons'].sum()
# Calculate the total number of codons inside hotspot regions (HR)
total_in_HR = disordered_codons['HR_codons'].sum()
# Calculate the total number of codons outside hotspot regions (HR)
total_out_HR = uniprot['CDS_aa'].sum() - disordered_codons['HR_codons'].sum(

print(f'Number of disorder codons out of HR: {disordered_out_HR}')
print(f'Number of disorder codons inside HR: {disordered_in_HR}')
print(f'Number of codons inside HR: {total_in_HR}')
print(f'Number of codons out of HR: {total_out_HR}')

proportion_in_HR = disordered_in_HR/total_in_HR
proportion_out_HR = disordered_out_HR/total_out_HR

# Fraction
print(f'Fraction of disorder codons inside HR: {round(proportion_in_HR,3)}')
print(f'Fraction of disorder codons out of HR: {round(proportion_out_HR,3)}')

# Percentage
print(f'Percentage of disorder codons inside HR: {round((proportion_in_HR*100),1)}')
print(f'Percentage of disorder codons out of HR: {round((proportion_out_HR*100),1)}')

from statsmodels.stats.proportion import proportions_ztest, proportion_confint

# Input data: number of disorder in each group and the total size of each group
count = [disordered_in_HR, disordered_out_HR]
nobs = [total_in_HR, total_out_HR]

# Z-test
stat, p_value = proportions_ztest(count, nobs)

# Confidence interval for the difference in proportions
confint = proportion_confint(count, nobs, alpha=0.05)

# Observed difference in proportions
prop_diff = (disordered_in_HR / total_in_HR) - (disordered_out_HR / total_out_HR)

# Confidence interval for the observed difference
ci_lower = prop_diff - (confint[0] - confint[1])[0]
ci_upper = prop_diff + (confint[0] - confint[1])[0]

# Print results and analysis
print(f'Z-statistic: {stat:.4f}')
print(f'P-value: {p_value:.4e}')
print(f'Difference in proportions: {prop_diff:.4f}')
print(f'95% CI for the difference in proportions: [{ci_lower:.4f}, {ci_upper:.4f}]')

