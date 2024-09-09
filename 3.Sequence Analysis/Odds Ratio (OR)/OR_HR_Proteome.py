# ## Packages
import re
import gzip
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

### 1. Data Preparation

# ### Hotspot Regions Sequences
# - DataFrame `df_HD` contains hotspot detection information with columns: 'Gene', 'Transcript', 'HR_coord', 'HR_len'
# - DataFrame `df_seq` contains sequences with columns: 'Transcript', 'AminoAcids' (full gene sequences)

# Merge df_HD and df_seq to obtain HR sequences for each gene
df_HR = pd.merge(df_HD, df_seq, on=['Transcript'], how='left')

# Function to extract the HR sequence based on coordinates
def extract_hr_seq(row):
    amino_acids = row['AminoAcids']
    htspt_coord = str(row['Htspt_coord'])  # Convert to string
    coords = htspt_coord.split(',')
    hr_seq_parts = []
    
    for coord in coords:
        start, end = map(int, coord.split('-'))
        start -= 1  # Adjust for 0-based indexing
        hr_seq_parts.append(amino_acids[start:end])
    
    hr_seq = ''.join(hr_seq_parts)
    return hr_seq

# Apply the function to extract HR sequences and add to a new column 'HR_seq'
df_HR['HR_seq'] = df_HR.apply(extract_hr_seq, axis=1)

# Drop the 'AminoAcids' column as it's no longer needed
df_HR.drop(columns=['AminoAcids'], inplace=True)

# Display the first few rows of the modified DataFrame
df_HR.head(5)

### 2. Amino Acid Counts

# #### Hotspot Region Amino Acid Counts
# - Count occurrences of each amino acid (AA) in each Hotspot Region (HR) sequence.
# - Then, sum the total occurrences of each AA in all HRs.

# List of amino acids
amino_acids = ['G', 'A', 'V', 'L', 'I', 'T', 'S', 'M', 'C', 'P', 'F', 'Y', 'W', 'H', 'K', 'R', 'D', 'E', 'N', 'Q', 'X']

# Function to count amino acid occurrences in each HR sequence
def count_amino_acids(hr_seq):
    counts = {aa: 0 for aa in amino_acids}
    sequences = hr_seq.split(',')
    for sequence in sequences:
        for aa in amino_acids:
            counts[aa] += sequence.count(aa)
    return counts

# Apply the function and create a DataFrame with AA counts
df_HR_AA = pd.concat([df_HR, df_HR['HR_seq'].apply(lambda x: pd.Series(count_amino_acids(x)))], axis=1)

# Sum all AAs and compare with 'HR_len'
df_HR_AA['Total'] = df_HR_AA[amino_acids].sum(axis=1)

# Sum across all HRs
HR_sum_counts = df_HR_AA[amino_acids].sum()
HR_sum_total = HR_sum_counts.sum()

# Create a summary DataFrame for HR amino acid counts
df_summary_HR_AA = pd.DataFrame(HR_sum_counts, columns=['counts']).reset_index()
df_summary_HR_AA = df_summary_HR_AA.rename(columns={'index': 'AA'})

# #### Proteome Amino Acid Counts
# Repeat the amino acid counting process for the full proteome sequences.

# Function to count amino acids in proteome sequences
df_proteome_AA = pd.concat([df_proteome, df_proteome['Sequence'].apply(lambda x: pd.Series(count_amino_acids(x)))], axis=1)

# Sum all AAs for the proteome
df_proteome_AA['Total'] = df_proteome_AA[amino_acids].sum(axis=1)
Prot_sum_counts = df_proteome_AA[amino_acids].sum()
Prot_sum_total = Prot_sum_counts.sum()

# Create a summary DataFrame for proteome amino acid counts
df_summary_proteome_AA = pd.DataFrame(Prot_sum_counts, columns=['counts']).reset_index()
df_summary_proteome_AA = df_summary_proteome_AA.rename(columns={'index': 'AA'})

### 3. Odds Calculation

# #### Hotspot Region Odds
df_summary_HR_AA['odds'] = df_summary_HR_AA['counts'] / (HR_sum_total - df_summary_HR_AA['counts'])

# #### Proteome Odds
df_summary_proteome_AA['odds'] = df_summary_proteome_AA['counts'] / (Prot_sum_total - df_summary_proteome_AA['counts'])

### 4. Odds Ratio and Confidence Interval Calculation

# Odds ratio (HR odds / Proteome odds)
or_hr_proteome = pd.merge(df_summary_HR_AA[['AA', 'odds']], df_summary_proteome_AA[['AA', 'odds']], on='AA', suffixes=('_HR', '_proteome'))
or_hr_proteome['odds_ratio'] = or_hr_proteome['odds_HR'] / or_hr_proteome['odds_proteome']

# Standard error of the log(OR)
se_log_or = np.sqrt((1 / or_hr_proteome['odds_HR']) + (1 / (HR_sum_total - or_hr_proteome['odds_HR'])) +
                    (1 / or_hr_proteome['odds_proteome']) + (1 / (Prot_sum_total - or_hr_proteome['odds_proteome'])))

# Z-value for 95% confidence intervals
z = 1.96

# Confidence intervals
or_hr_proteome['ci_lower'] = np.exp(np.log(or_hr_proteome['odds_ratio']) - z * se_log_or)
or_hr_proteome['ci_upper'] = np.exp(np.log(or_hr_proteome['odds_ratio']) + z * se_log_or)

# Display results
or_hr_proteome[['AA', 'odds_ratio', 'ci_lower', 'ci_upper']]

### 5. Plotting Odds Ratios

# Sort by odds_ratio in descending order
or_hr_proteome = or_hr_proteome.sort_values(by='odds_ratio', ascending=False)

# Set up the plot
fig, ax = plt.subplots(figsize=(10, 6))

# Colors for bars (different colors for odds_ratio >= 1 and < 1)
colors = ['#4d82bc' if val >= 1 else '#c4dafa' for val in or_hr_proteome['odds_ratio']]

# Plot the bars
for idx, (label, odds_ratio) in enumerate(zip(or_hr_proteome['AA'], or_hr_proteome['odds_ratio'])):
    if odds_ratio >= 1:
        ax.bar(label, odds_ratio - 1, bottom=1, color=colors[idx])
    else:
        ax.bar(label, 1 - odds_ratio, bottom=odds_ratio, color=colors[idx])

# Add labels and confidence intervals
for idx, (label, odds_ratio) in enumerate(zip(or_hr_proteome['AA'], or_hr_proteome['odds_ratio'])):
    if odds_ratio >= 1:
        ax.text(label, 1 + (odds_ratio - 1) + 0.1, f'{odds_ratio:.2f}', ha='center', va='bottom')
    else:
        ax.text(label, odds_ratio - 0.1, f'{odds_ratio:.2f}', ha='center', va='top')

# Add confidence intervals
ax.errorbar(or_hr_proteome['AA'], or_hr_proteome['odds_ratio'], 
            yerr=[or_hr_proteome['odds_ratio'] - or_hr_proteome['ci_lower'], or_hr_proteome['ci_upper'] - or_hr_proteome['odds_ratio']], 
            fmt='none', ecolor='black', capsize=3, capthick=0.5, linewidth=0.5)

# Add a dashed line at y=1
ax.axhline(y=1, color='red', linestyle='--', linewidth=1.5)

# Set labels and title
ax.set_xlabel('Amino Acids', labelpad=20)
ax.set_ylabel('Odds Ratios', labelpad=20)

# Adjust the y-axis range for better visualization
max_y = or_hr_proteome['odds_ratio'].max() + 0.5
ax.set_ylim(-0.25, max_y)

# Show the plot
plt.show()
