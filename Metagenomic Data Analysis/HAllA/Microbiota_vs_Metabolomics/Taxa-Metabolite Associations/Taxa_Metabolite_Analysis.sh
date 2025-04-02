
#!/bin/bash

# Activate the Python environment (if needed)
# source activate my_env

# Run the Python analysis script
python3 <<EOF

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages

# Load the input files
associations_file = 'Associations.xlsx'
metabolites_file = 'Metabolitesv2.xlsx'
tsv_file = 'X_original.tsv'

# Load data from the Excel files
associations_df = pd.read_excel(associations_file, sheet_name='Sheet1')
metabolites_df = pd.read_excel(metabolites_file, sheet_name='Sheet1')
tsv_data = pd.read_csv(tsv_file, sep='\t')

# Analysis 1: Generate Taxa vs Metabolite Associations Bar Plot
taxa_metabolite_counts = associations_df.groupby('Taxa')['Metabolite'].nunique().sort_values()

plt.figure(figsize=(15, 8))
plt.bar(taxa_metabolite_counts.index, taxa_metabolite_counts.values)
plt.xticks(rotation=90, fontsize=8)
plt.xlabel("Individual Taxa (sorted by number of associations)")
plt.ylabel("Total Number of Metabolite Associations")
plt.title("Taxa Stratified by Metabolite Associations")
plt.tight_layout()

# Save PDF
taxa_pdf_file = 'Taxa_Metabolite_Associations.pdf'
with PdfPages(taxa_pdf_file) as pdf:
    pdf.savefig()
    plt.close()

# Analysis 2: Generate Dot-Line Plot
plt.figure(figsize=(15, 8))
plt.plot(taxa_metabolite_counts.index, taxa_metabolite_counts.values, marker='o', linestyle='-', color='b')
plt.xticks(rotation=90, fontsize=8)
plt.xlabel("Individual Taxa (sorted by number of associations)")
plt.ylabel("Total Number of Metabolite Associations")
plt.title("Dot-Line Plot: Taxa Stratified by Metabolite Associations")
plt.tight_layout()

# Save PDF
dotline_pdf_file = 'Taxa_Metabolite_DotLine_Plot.pdf'
with PdfPages(dotline_pdf_file) as pdf:
    pdf.savefig()
    plt.close()

# Analysis 3: Metabolite Associations Analysis
limited_metabolite_taxa_counts = associations_df.groupby('Metabolite')['Taxa'].nunique()
sorted_metabolite_taxa_counts = limited_metabolite_taxa_counts.sort_values()

# Dot-Line Plot for Metabolites
plt.figure(figsize=(15, 8))
plt.plot(sorted_metabolite_taxa_counts.index, sorted_metabolite_taxa_counts.values, marker='o', linestyle='-', color='g')
plt.xticks(rotation=90, fontsize=8)
plt.xlabel("Individual Metabolites (sorted by number of taxa associations)")
plt.ylabel("Total Number of Taxa Associations")
plt.title("Dot-Line Plot: Metabolites Stratified by Number of Taxa Associations")
plt.tight_layout()

# Save PDF
metabolite_dotline_pdf_file = 'Metabolite_Taxa_DotLine_Plot.pdf'
with PdfPages(metabolite_dotline_pdf_file) as pdf:
    pdf.savefig()
    plt.close()

# Save Excel file with metabolite-taxa associations
metabolite_taxa_file = 'Metabolite_Taxa_Distribution.xlsx'
metabolite_taxa_df = associations_df[['Taxa', 'Metabolite']].drop_duplicates()
metabolite_taxa_df.to_excel(metabolite_taxa_file, index=False, sheet_name='Distribution_Bins')

# Analysis 4: Generate Heatmap of Taxa Presence Across Vivaria
limited_taxa_list = metabolite_taxa_df['Taxa'].unique()
filtered_taxa_abundance = tsv_data[tsv_data['MicrobiotaName'].isin(limited_taxa_list)]

# Convert positive values to 1 (presence) and zero values to 0 (absence)
binary_taxa_abundance = filtered_taxa_abundance.set_index('MicrobiotaName').gt(0).astype(int)

# Create Heatmap
plt.figure(figsize=(20, 12))
sns.heatmap(binary_taxa_abundance, cmap="YlGnBu", cbar_kws={'label': 'Presence (1) / Absence (0)'}, linewidths=.5)
plt.title("Corrected Presence-Absence Heatmap of Limited Taxa Across Vivaria")
plt.xlabel("Vivaria Samples")
plt.ylabel("Taxa")
plt.tight_layout()

# Save PDF
heatmap_pdf_file = 'Limited_Taxa_Presence_Heatmap.pdf'
with PdfPages(heatmap_pdf_file) as pdf:
    pdf.savefig()
    plt.close()

# Print Completion Message
print("Analysis completed. Files generated:")
print(f"1. {taxa_pdf_file}")
print(f"2. {dotline_pdf_file}")
print(f"3. {metabolite_dotline_pdf_file}")
print(f"4. {metabolite_taxa_file} (Excel)")
print(f"5. {heatmap_pdf_file}")

EOF
