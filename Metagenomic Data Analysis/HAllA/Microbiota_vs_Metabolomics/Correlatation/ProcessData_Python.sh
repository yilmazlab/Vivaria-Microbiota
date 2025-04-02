
#!/bin/bash

# Activate Python environment if necessary
# source activate my_environment

# Execute Python script
python <<EOF
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load the data
metabolite_list = pd.read_excel("MetaboliteList.xlsx", header=None)[0].tolist()
associations_data = pd.read_excel("Associations.xlsx")
x_original_data = pd.read_csv("X_original.tsv", sep='\t')

# Initialize a dictionary for concatenated abundance data
concatenated_abundance_data = {}

# Iterate over each metabolite
for metabolite in metabolite_list:
    associated_bacteria = associations_data[associations_data['Metabolite'] == metabolite]['Taxa']
    
    abundance_values = []
    for bacterium in associated_bacteria:
        if bacterium in x_original_data['MicrobiotaName'].values:
            abundance = x_original_data.loc[x_original_data['MicrobiotaName'] == bacterium].iloc[:, 1:]
            abundance_values.append(abundance.squeeze())

    concatenated_abundance_data[metabolite] = pd.concat(abundance_values, ignore_index=True)

# Convert to DataFrame
concatenated_abundance_df = pd.DataFrame(concatenated_abundance_data)

# Save to Excel
concatenated_abundance_df.to_excel("ConcatenatedRelativeAbundanceData.xlsx", index=False)

# Generate dot plot
plt.figure(figsize=(20, 10))
for idx, column in enumerate(concatenated_abundance_df.columns):
    plt.scatter([idx] * len(concatenated_abundance_df[column]), np.log10(concatenated_abundance_df[column] + 1e-10), color='blue', alpha=0.6)
plt.xticks(range(len(concatenated_abundance_df.columns)), concatenated_abundance_df.columns, rotation=90)
plt.ylabel('Log10(Relative Abundance)')
plt.title('Dot Plot of Concatenated Relative Abundance')
plt.tight_layout()
plt.savefig("ConcatenatedRelativeAbundanceDotPlot.pdf", format='pdf')
EOF
